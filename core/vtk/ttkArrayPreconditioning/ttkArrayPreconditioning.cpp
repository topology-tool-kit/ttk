#include <OrderDisambiguation.h>
#include <Triangulation.h>
#include <ttkArrayPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <regex>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkArrayPreconditioning);

ttkArrayPreconditioning::ttkArrayPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkArrayPreconditioning::Modified);
}

int ttkArrayPreconditioning::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkArrayPreconditioning::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkArrayPreconditioning::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  ttk::Timer tm{};

  int keepGoing = checkEmptyMPIInput<vtkDataSet>(input);
  if(keepGoing < 2) {
    return keepGoing;
  }
#ifdef TTK_ENABLE_MPI
  const auto triangulation{this->GetTriangulation(input)};
#endif
  output->ShallowCopy(input);

  auto pointData = input->GetPointData();
  size_t nVertices = input->GetNumberOfPoints();

  std::vector<vtkDataArray *> scalarArrays{};

  if(SelectFieldsWithRegexp) {
    // select all input point data arrays whose name is matching the regexp
    const auto n = pointData->GetNumberOfArrays();
    for(int i = 0; i < n; ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && std::regex_match(array->GetName(), std::regex(RegexpString))) {
        scalarArrays.emplace_back(array);
      }
    }
  } else {
    // get all selected input point data arrays
    for(int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && ArraySelection->ArrayIsEnabled(array->GetName())) {
        scalarArrays.emplace_back(array);
      }
    }
  }

#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    // add the order array for every scalar array, except the ghostcells, the
    // rankarray and the global ids
    for(auto scalarArray : scalarArrays) {
      int status = 0;
      std::string arrayName = std::string(scalarArray->GetName());
      if(arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"
         && arrayName != "RankArray") {
        this->printMsg("Arrayname: " + arrayName);
        vtkNew<ttkSimplexIdTypeArray> orderArray{};
        orderArray->SetName(
          ttkArrayPreconditioning::GetOrderArrayName(scalarArray).data());
        orderArray->SetNumberOfComponents(1);
        orderArray->SetNumberOfTuples(nVertices);

        this->printMsg(std::to_string(scalarArray->GetDataType()));
        ttkTypeMacroA(
          scalarArray->GetDataType(),
          (status = processScalarArray(
             ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
             ttkUtils::GetPointer<T0>(scalarArray),
             [triangulation](const ttk::SimplexId a) {
               return triangulation->getVertexGlobalId(a);
             },
             [triangulation](const ttk::SimplexId a) {
               return triangulation->getVertexRank(a);
             },
             [triangulation](const ttk::SimplexId a) {
               return triangulation->getVertexLocalId(a);
             },
             nVertices, BurstSize, triangulation->getNeighborRanks())));

        // On error cancel filter execution
        if(status != 1)
          return 0;
        output->GetPointData()->AddArray(orderArray);
      }
    }
    this->printMsg("Preconditioned selected scalar arrays", 1.0,
                   tm.getElapsedTime(), this->threadNumber_);
    return 1;
  } else {
    this->printMsg("Necessary arrays are present,  TTK is built with MPI "
                   "support, but not run with mpirun. Running sequentially.");
  }
#else
  this->printMsg("Necessary arrays are present, but TTK is not built with "
                 "MPI support, running sequentially.");
#endif

  for(auto scalarArray : scalarArrays) {
    vtkNew<ttkSimplexIdTypeArray> orderArray{};
    orderArray->SetName(
      ttkArrayPreconditioning::GetOrderArrayName(scalarArray).data());
    orderArray->SetNumberOfComponents(1);
    orderArray->SetNumberOfTuples(nVertices);

    switch(scalarArray->GetDataType()) {
      vtkTemplateMacro(ttk::sortVertices(
        nVertices, static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
        static_cast<int *>(nullptr),
        static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray)),
        this->threadNumber_));
    }

    output->GetPointData()->AddArray(orderArray);
    this->printMsg("Generated order array for scalar array `"
                   + std::string{scalarArray->GetName()} + "'");
  }

  this->printMsg("Preconditioned selected scalar arrays", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
