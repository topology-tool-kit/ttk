#include <ttkDepthImageBasedGeometryApproximation.h>

#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDepthImageBasedGeometryApproximation)

  int ttkDepthImageBasedGeometryApproximation::RequestData(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
  Memory mem;
  Timer t;

  // Print status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkDepthImageBasedGeometryApproximation] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Set Wrapper
  depthImageBasedGeometryApproximation_.setWrapper(this);

  // Prepare input and output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto inputMBD = vtkMultiBlockDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto outputMBD = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Process each depth image individually
  size_t nBlocks = inputMBD->GetNumberOfBlocks();
  for(size_t i = 0; i < nBlocks; i++) {
    // Get vtkImageData
    auto inputImage = vtkImageData::SafeDownCast(inputMBD->GetBlock(i));

    // Get input paramters
    auto inPointData = inputImage->GetPointData();
    auto depthValues = inputImage->GetPointData()->GetAbstractArray(
      this->GetDepthScalarField().data());

    auto inFieldData = inputImage->GetFieldData();
    auto camHeight = inFieldData->GetAbstractArray("CamHeight");
    auto camPosition = inFieldData->GetAbstractArray("CamPosition");
    auto camDirection = inFieldData->GetAbstractArray("CamDirection");
    auto camUp = inFieldData->GetAbstractArray("CamUp");
    auto camNearFar = inFieldData->GetAbstractArray("CamNearFar");
    auto camRes = inFieldData->GetAbstractArray("CamRes");

    // Check if all parameters are present
    if(depthValues == nullptr || camHeight == nullptr || camPosition == nullptr
       || camDirection == nullptr || camUp == nullptr || camNearFar == nullptr
       || camRes == nullptr) {
      stringstream msg;
      msg << "[ttkDepthImageBasedGeometryApproximation] ERROR: Input depth "
             "image does not have one or more of the required fields (see "
             "Cinema Spec D - Data Product Specification)"
          << endl;
      dMsg(cout, msg.str(), fatalMsg);
      return 0;
    }

    // Vectors to hold raw approximated geometry
    vector<size_t> indicies;
    vector<tuple<double, double, double>> vertices;
    vector<tuple<int, int, int>> triangles;
    vector<double> triangleDistortions;

    // Approximate geometry
    switch(depthValues->GetDataType()) {
      vtkTemplateMacro({
        depthImageBasedGeometryApproximation_.execute<VTK_TT>(
          // Input
          (VTK_TT *)depthValues->GetVoidPointer(0),
          (double *)camPosition->GetVoidPointer(0),
          (double *)camDirection->GetVoidPointer(0),
          (double *)camUp->GetVoidPointer(0),
          (double *)camRes->GetVoidPointer(0),
          (double *)camNearFar->GetVoidPointer(0),
          (double *)camHeight->GetVoidPointer(0), this->GetSubsampling(),

          // Output
          indicies, vertices, triangles, triangleDistortions);
      });
    }

    // Represent approximated geometry via VTK
    auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create points
    {
      size_t n = vertices.size();

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(n);

      auto pointCoords = (float *)points->GetVoidPointer(0);
      size_t j = 0;
      for(auto &x : vertices) {
        pointCoords[j++] = get<0>(x);
        pointCoords[j++] = get<1>(x);
        pointCoords[j++] = get<2>(x);
      }

      mesh->SetPoints(points);
    }

    // Copy Point Data
    {
      size_t n = inPointData->GetNumberOfArrays();

      auto outPointData = mesh->GetPointData();
      size_t m = indicies.size();

      for(size_t j = 0; j < n; j++) {
        auto inArray = inPointData->GetArray(j);

        auto outArray = vtkDataArray::CreateDataArray(inArray->GetDataType());
        outArray->SetName(inArray->GetName());
        outArray->SetNumberOfTuples(m);

        for(size_t k = 0; k < m; k++) {
          outArray->SetTuple(k, inArray->GetTuple(indicies[k]));
        }

        outPointData->AddArray(outArray);
      }
    }

    // Create cells
    {
      size_t n = triangles.size();

      auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
      cells->SetNumberOfValues(4 * n);
      auto cellIds = (vtkIdType *)cells->GetVoidPointer(0);

      auto triangleDistortionsScalars = vtkSmartPointer<vtkDoubleArray>::New();
      triangleDistortionsScalars->SetNumberOfValues(n);
      triangleDistortionsScalars->SetNumberOfComponents(1);
      triangleDistortionsScalars->SetName("Distortion");
      double *triangleDistortionsScalarsData
        = (double *)triangleDistortionsScalars->GetVoidPointer(0);

      size_t q = 0;
      for(size_t j = 0; j < triangles.size(); j++) {
        cellIds[q++] = 3;
        cellIds[q++] = (vtkIdType)get<0>(triangles[j]);
        cellIds[q++] = (vtkIdType)get<1>(triangles[j]);
        cellIds[q++] = (vtkIdType)get<2>(triangles[j]);

        triangleDistortionsScalarsData[j] = triangleDistortions[j];
      }

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetCells(n, cells);
      mesh->SetCells(VTK_TRIANGLE, cellArray);

      mesh->GetCellData()->AddArray(triangleDistortionsScalars);
    }

    outputMBD->SetBlock(i, mesh);
    this->updateProgress(((float)i) / ((float)(nBlocks - 1)));
  }

  // Print status
  {
    stringstream msg;
    msg << "[ttkDepthImageBasedGeometryApproximation] "
           "--------------------------------------"
        << endl
        << "[ttkDepthImageBasedGeometryApproximation] " << nBlocks
        << " Images processed" << endl
        << "[ttkDepthImageBasedGeometryApproximation]   Time: "
        << t.getElapsedTime() << " s" << endl
        << "[ttkDepthImageBasedGeometryApproximation] Memory: "
        << mem.getElapsedUsage() << " MB" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
