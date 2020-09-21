#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkPersistenceDiagram.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagram);

ttkPersistenceDiagram::ttkPersistenceDiagram() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistenceDiagram::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagram(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<ttk::PersistencePair> &diagram,
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<ttkSimplexIdTypeArray> vertexIdentifierScalars{};
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkNew<vtkIntArray> nodeTypeScalars{};
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkNew<ttkSimplexIdTypeArray> pairIdentifierScalars{};
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> extremumIndexScalars{};
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkNew<vtkFloatArray> coordsScalars{};
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];
    vtkIdType oldIds[2];

    double maxPersistenceValue = std::numeric_limits<double>::min();
    oldIds[0] = 0;
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const double persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      std::array<double, 3> p{0, 0, 0};
      const auto a = std::get<0>(diagram[i]);
      const auto na = static_cast<ttk::SimplexId>(std::get<1>(diagram[i]));
      const auto b = std::get<2>(diagram[i]);
      const auto nb = static_cast<ttk::SimplexId>(std::get<3>(diagram[i]));

      nodeTypeScalars->InsertTuple1(2 * i, na);
      nodeTypeScalars->InsertTuple1(2 * i + 1, nb);

      vertexIdentifierScalars->InsertTuple1(2 * i, a);
      vertexIdentifierScalars->InsertTuple1(2 * i + 1, b);

      std::array<float, 3> coords{};
      triangulation->getVertexPoint(a, coords[0], coords[1], coords[2]);
      coordsScalars->InsertTuple3(2 * i, coords[0], coords[1], coords[2]);

      triangulation->getVertexPoint(b, coords[0], coords[1], coords[2]);
      coordsScalars->InsertTuple3(2 * i + 1, coords[0], coords[1], coords[2]);

      p[0] = inputScalars->GetTuple1(a);
      p[1] = inputScalars->GetTuple1(a);
      ids[0] = points->InsertNextPoint(p.data());

      p[0] = inputScalars->GetTuple1(a);
      p[1] = inputScalars->GetTuple1(b);
      ids[1] = points->InsertNextPoint(p.data());

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
    oldIds[1] = ids[0];

    // add diag
    persistenceDiagram->InsertNextCell(VTK_LINE, 2, oldIds);
    pairIdentifierScalars->InsertTuple1(diagramSize, -1);
    extremumIndexScalars->InsertTuple1(diagramSize, -1);
    persistenceScalars->InsertTuple1(diagramSize, 2 * maxPersistenceValue);
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);

  return 0;
}

template <class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagramInsideDomain(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<ttk::PersistencePair> &diagram,
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<ttkSimplexIdTypeArray> vertexIdentifierScalars{};
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkNew<vtkIntArray> nodeTypeScalars{};
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkNew<ttkSimplexIdTypeArray> pairIdentifierScalars{};
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> extremumIndexScalars{};
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkDataArray *birthScalars = inputScalars->NewInstance();
  birthScalars->SetNumberOfComponents(1);
  birthScalars->SetName("Birth");

  vtkDataArray *deathScalars = inputScalars->NewInstance();
  deathScalars->SetNumberOfComponents(1);
  deathScalars->SetName("Death");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];

    double maxPersistenceValue = std::numeric_limits<double>::min();
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const double persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      std::array<float, 3> p{};
      const auto a = std::get<0>(diagram[i]);
      const auto na = static_cast<ttk::SimplexId>(std::get<1>(diagram[i]));
      const auto b = std::get<2>(diagram[i]);
      const auto nb = static_cast<ttk::SimplexId>(std::get<3>(diagram[i]));
      const double sa = inputScalars->GetTuple1(a);
      const double sb = inputScalars->GetTuple1(b);

      nodeTypeScalars->InsertTuple1(2 * i, na);
      nodeTypeScalars->InsertTuple1(2 * i + 1, nb);
      vertexIdentifierScalars->InsertTuple1(2 * i, a);
      vertexIdentifierScalars->InsertTuple1(2 * i + 1, b);
      birthScalars->InsertTuple1(2 * i, sa);
      birthScalars->InsertTuple1(2 * i + 1, sa);
      deathScalars->InsertTuple1(2 * i, sa);
      deathScalars->InsertTuple1(2 * i + 1, sb);

      triangulation->getVertexPoint(a, p[0], p[1], p[2]);
      ids[0] = points->InsertNextPoint(p.data());

      triangulation->getVertexPoint(b, p[0], p[1], p[2]);
      ids[1] = points->InsertNextPoint(p.data());

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(birthScalars);
  persistenceDiagram->GetPointData()->AddArray(deathScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);
  return 0;
}

void ttkPersistenceDiagram::Modified() {
  computeDiagram_ = true;
  ttkAlgorithm::Modified();
}

int ttkPersistenceDiagram::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputCTPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars) {
    this->printErr("Wrong input scalars");
    return 0;
  }
#endif

  vtkDataArray *offsetField
    = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!offsetField) {
    this->printErr("Wrong input offsets");
    return 0;
  }
  if(offsetField->GetDataType() != VTK_INT
     and offsetField->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported");
    return 0;
  }
#endif

  if(this->GetMTime() < inputScalars->GetMTime())
    computeDiagram_ = true;

  int status = 0;
  if(computeDiagram_) {
    CTDiagram_.clear();
    ttkVtkTemplateMacro(
      inputScalars->GetDataType(), triangulation->getType(),
      status = this->execute(
        CTDiagram_,
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
        static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsetField)),
        static_cast<TTK_TT *>(triangulation->getData())));

    // something wrong in baseCode
    if(status) {
      std::stringstream msg;
      msg << "PersistenceDiagram::execute() error code : " << status;
      this->printErr(msg.str());
      return 0;
    }
  }

  if(ShowInsideDomain) {
    ttkTemplateMacro(
      triangulation->getType(),
      setPersistenceDiagramInsideDomain(
        outputCTPersistenceDiagram, ttk::ftm::TreeType::Contour, CTDiagram_,
        inputScalars, static_cast<TTK_TT *>(triangulation->getData())));
  } else {
    ttkTemplateMacro(
      triangulation->getType(),
      setPersistenceDiagram(
        outputCTPersistenceDiagram, ttk::ftm::TreeType::Contour, CTDiagram_,
        inputScalars, static_cast<TTK_TT *>(triangulation->getData())));
  }

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  computeDiagram_ = false;

  return 1;
}
