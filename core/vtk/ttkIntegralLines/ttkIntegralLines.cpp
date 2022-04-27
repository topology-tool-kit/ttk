#include <ttkIntegralLines.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <array>

vtkStandardNewMacro(ttkIntegralLines);

ttkIntegralLines::ttkIntegralLines() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkIntegralLines::~ttkIntegralLines() = default;

int ttkIntegralLines::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkIntegralLines::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

int ttkIntegralLines::getTrajectories(
  vtkDataSet *input,
  ttk::Triangulation *triangulation,
  std::vector<std::vector<ttk::SimplexId>> &trajectories,
  vtkUnstructuredGrid *output) {

  if(input == nullptr || output == nullptr
     || input->GetPointData() == nullptr) {
    this->printErr("Null pointers in getTrajectories parameters");
    return 0;
  }

  vtkNew<vtkUnstructuredGrid> ug{};
  vtkNew<vtkPoints> pts{};

  vtkNew<vtkFloatArray> dist{};
  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");

  const auto numberOfArrays = input->GetPointData()->GetNumberOfArrays();

  std::vector<vtkDataArray *> scalarArrays{};
  scalarArrays.reserve(numberOfArrays);
  for(int k = 0; k < numberOfArrays; ++k) {
    const auto a = input->GetPointData()->GetArray(k);
    if(a->GetNumberOfComponents() == 1) {
      // only keep scalar arrays
      scalarArrays.push_back(a);
    }
  }

  std::vector<vtkSmartPointer<vtkDataArray>> inputScalars(scalarArrays.size());
  for(size_t k = 0; k < scalarArrays.size(); ++k) {
    inputScalars[k]
      = vtkSmartPointer<vtkDataArray>::Take(scalarArrays[k]->NewInstance());
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }

  std::array<float, 3> p0{}, p1{};
  std::array<vtkIdType, 2> ids{};
  for(size_t i = 0; i < trajectories.size(); ++i) {
    if(!trajectories[i].empty()) {
      auto vertex = trajectories[i][0];
      // init
      triangulation->getVertexPoint(vertex, p0[0], p0[1], p0[2]);
      ids[0] = pts->InsertNextPoint(p0.data());
      // distanceScalars
      float distanceFromSeed{};
      dist->InsertNextTuple1(distanceFromSeed);
      // inputScalars
      for(size_t k = 0; k < scalarArrays.size(); ++k) {
        inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));
      }

      for(size_t j = 1; j < trajectories[i].size(); ++j) {
        vertex = trajectories[i][j];
        triangulation->getVertexPoint(vertex, p1[0], p1[1], p1[2]);
        ids[1] = pts->InsertNextPoint(p1.data());
        // distanceScalars
        distanceFromSeed += ttk::Geometry::distance(p0.data(), p1.data(), 3);
        dist->InsertNextTuple1(distanceFromSeed);
        // inputScalars
        for(size_t k = 0; k < scalarArrays.size(); ++k) {
          inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));
        }

        ug->InsertNextCell(VTK_LINE, 2, ids.data());

        // iteration
        ids[0] = ids[1];
        p0 = p1;
      }
    }
  }
  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  for(size_t k = 0; k < scalarArrays.size(); ++k) {
    ug->GetPointData()->AddArray(inputScalars[k]);
  }

  output->ShallowCopy(ug);

  return 1;
}

int ttkIntegralLines::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0], 0);
  vtkPointSet *seeds = vtkPointSet::GetData(inputVector[1], 0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(domain);
  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, domain);

  vtkDataArray *inputOffsets
    = this->GetOrderArray(domain, 0, 1, ForceInputOffsetScalarField);

  std::vector<ttk::SimplexId> idSpareStorage{};
  auto *inputIdentifiers = this->GetIdentifierArrayPtr(
    ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
    idSpareStorage);

  const auto numberOfPointsInDomain = domain->GetNumberOfPoints();
  const auto numberOfPointsInSeeds = seeds->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  // triangulation problem
  if(!triangulation) {
    this->printErr("wrong triangulation.");
    return -1;
  }
  // field problem
  if(!inputScalars) {
    this->printErr("wrong scalar field.");
    return -1;
  }
  // field problem
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }
  // field problem
  if(!inputIdentifiers) {
    this->printErr("wrong identifiers.");
    return -1;
  }
  // no points.
  if(numberOfPointsInDomain <= 0) {
    this->printErr("domain has no points.");
    return -1;
  }
  // no points.
  if(numberOfPointsInSeeds <= 0) {
    this->printErr("seeds have no points.");
    return -1;
  }
#endif

  std::vector<std::vector<ttk::SimplexId>> trajectories{};

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSeedNumber(numberOfPointsInSeeds);
  this->setDirection(Direction);
  this->setInputScalarField(inputScalars->GetVoidPointer(0));
  this->setInputOffsets(ttkUtils::GetPointer<ttk::SimplexId>(inputOffsets));
  this->setVertexIdentifierScalarField(inputIdentifiers);
  this->setOutputTrajectories(&trajectories);

  this->preconditionTriangulation(triangulation);

  int status = 0;
  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (status = this->execute<VTK_TT, TTK_TT>(
                         static_cast<TTK_TT *>(triangulation->getData()))));
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(status != 0) {
    this->printErr("IntegralLines.execute() error code : "
                   + std::to_string(status));
    return 0;
  }
#endif

  // make the vtk trajectories
  status = getTrajectories(domain, triangulation, trajectories, output);
  return status;
}
