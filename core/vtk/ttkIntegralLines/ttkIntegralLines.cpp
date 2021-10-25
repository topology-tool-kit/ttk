#include <ttkIntegralLines.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIntegralLines)

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

int ttkIntegralLines::getTrajectories(vtkDataSet *input,
                                      ttk::Triangulation *triangulation,
                                      vector<vector<SimplexId>> &trajectories,
                                      vtkUnstructuredGrid *output) {
  vtkSmartPointer<vtkUnstructuredGrid> ug
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> dist = vtkSmartPointer<vtkFloatArray>::New();
  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");

  // here, copy the original scalars
  int numberOfArrays = input->GetPointData()->GetNumberOfArrays();

  vector<vtkDataArray *> scalarArrays;
  for(int k = 0; k < numberOfArrays; ++k) {
    auto a = input->GetPointData()->GetArray(k);

    if(a->GetNumberOfComponents() == 1)
      scalarArrays.push_back(a);
  }
  // not efficient, implicit conversion to double
  vector<vtkSmartPointer<vtkDoubleArray>> inputScalars(scalarArrays.size());
  for(unsigned int k = 0; k < scalarArrays.size(); ++k) {
    inputScalars[k] = vtkSmartPointer<vtkDoubleArray>::New();
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }

  float p0[3];
  float p1[3];
  vtkIdType ids[2];
  for(SimplexId i = 0; i < (SimplexId)trajectories.size(); ++i) {
    if(trajectories[i].size()) {
      SimplexId vertex = trajectories[i][0];
      // init
      triangulation->getVertexPoint(vertex, p0[0], p0[1], p0[2]);
      ids[0] = pts->InsertNextPoint(p0);
      // distanceScalars
      float distanceFromSeed{};
      dist->InsertNextTuple1(distanceFromSeed);
      // inputScalars
      for(unsigned int k = 0; k < scalarArrays.size(); ++k)
        inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

      for(SimplexId j = 1; j < (SimplexId)trajectories[i].size(); ++j) {
        vertex = trajectories[i][j];
        triangulation->getVertexPoint(vertex, p1[0], p1[1], p1[2]);
        ids[1] = pts->InsertNextPoint(p1);
        // distanceScalars
        distanceFromSeed += Geometry::distance(p0, p1, 3);
        dist->InsertNextTuple1(distanceFromSeed);
        // inputScalars
        for(unsigned int k = 0; k < scalarArrays.size(); ++k)
          inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

        ug->InsertNextCell(VTK_LINE, 2, ids);

        // iteration
        ids[0] = ids[1];
        p0[0] = p1[0];
        p0[1] = p1[1];
        p0[2] = p1[2];
      }
    }
  }
  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  for(unsigned int k = 0; k < scalarArrays.size(); ++k)
    ug->GetPointData()->AddArray(inputScalars[k]);

  output->ShallowCopy(ug);

  return 0;
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

  std::vector<SimplexId> idSpareStorage{};
  auto *inputIdentifiers = this->GetIdentifierArrayPtr(
    ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
    idSpareStorage);

  const SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
  const SimplexId numberOfPointsInSeeds = seeds->GetNumberOfPoints();

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

  vector<vector<SimplexId>> trajectories;

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSeedNumber(numberOfPointsInSeeds);
  this->setDirection(Direction);
  this->setInputScalarField(inputScalars->GetVoidPointer(0));
  this->setInputOffsets(
    static_cast<SimplexId *>(inputOffsets->GetVoidPointer(0)));

  this->setVertexIdentifierScalarField(inputIdentifiers);
  this->setOutputTrajectories(&trajectories);

  this->preconditionTriangulation(triangulation);

  int status = 0;
  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (status = this->execute<VTK_TT, TTK_TT>(
                         static_cast<TTK_TT *>(triangulation->getData()))));
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(status) {
    std::stringstream msg;
    msg << "IntegralLines.execute() error code : " << status;
    this->printErr(msg.str());
    return -1;
  }
#endif

  // make the vtk trajectories
  getTrajectories(domain, triangulation, trajectories, output);

  return (int)(status == 0);
}
