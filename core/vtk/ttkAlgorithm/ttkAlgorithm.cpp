#include <ttkAlgorithm.h>

#include <vtkDataObject.h> // For output port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkSmartPointer.h>

#include <vtkInformationVector.h>

#include <vtkCompositeDataPipeline.h>

#include <vtkFieldData.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// Pass input type information key
#include <vtkInformationKey.h>
vtkInformationKeyMacro(ttkAlgorithm, SAME_DATA_TYPE_AS_INPUT_PORT, Integer);

vtkStandardNewMacro(ttkAlgorithm);

ttkAlgorithm::ttkAlgorithm() {
}
ttkAlgorithm::~ttkAlgorithm() {
}

// init static triangulation registry
std::unordered_map<void *, std::pair<ttk::Triangulation *, vtkMTimeType>>
  ttkAlgorithm::DataSetToTriangulationMap;

ttk::Triangulation *ttkAlgorithm::GetTriangulation(vtkDataSet *dataSet) {

  this->PrintMsg("Requesting triangulation for '"
                   + std::string(dataSet->GetClassName()) + "'",
                 ttk::debug::Priority::DETAIL);

  switch(dataSet->GetDataObjectType()) {
    // =====================================================================
    case VTK_UNSTRUCTURED_GRID: {
      auto dataSetAsUG = (vtkUnstructuredGrid *)dataSet;

      // check if triangulation already exists
      auto it = ttkAlgorithm::DataSetToTriangulationMap.find(
        (void *)dataSetAsUG->GetCells());
      if(it != ttkAlgorithm::DataSetToTriangulationMap.end()) {
        if(it->second.second == dataSetAsUG->GetCells()->GetMTime()) {
          this->PrintMsg("Returning already initilized triangulation",
                         ttk::debug::Priority::DETAIL);
          return it->second.first;
        } else {
          this->PrintMsg("Chached triangulation no longer valid",
                         ttk::debug::Priority::DETAIL);
        }
      }

      this->PrintMsg(
        "Initializing triangulation", 0, ttk::debug::Priority::DETAIL);
      auto newTriangulation = new ttk::Triangulation();
      auto cells = dataSetAsUG->GetCells();

      // init points
      if(dataSet->GetNumberOfPoints() > 0) {
        auto points = dataSetAsUG->GetPoints();
        auto pointDataType = points->GetDataType();
        if(pointDataType != VTK_FLOAT && pointDataType != VTK_DOUBLE) {
          this->PrintErr("Unable to initialize 'ttk::Triangulation' for point "
                         "precision other than 'float' or 'double'.");
          delete newTriangulation;
          return nullptr;
        }

        newTriangulation->setInputPoints(dataSet->GetNumberOfPoints(),
                                         points->GetVoidPointer(0),
                                         pointDataType == VTK_DOUBLE);
      }

      // init cells
      if(dataSet->GetNumberOfCells() > 0) {
        newTriangulation->setInputCells(
          dataSet->GetNumberOfCells(), cells->GetPointer());
      }

      ttkAlgorithm::DataSetToTriangulationMap.insert(
        {(void *)cells, {newTriangulation, cells->GetMTime()}});

      this->PrintMsg("Initializing Triangulation", 1,
                     ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);
      return newTriangulation;
    }

    // =====================================================================
    case VTK_POLY_DATA: {
      auto dataSetAsPD = (vtkPolyData *)dataSet;

      // check if triangulation already exists
      auto it = ttkAlgorithm::DataSetToTriangulationMap.end();
      int triangulationStatus
        = 0; // 0: no triangulation chached; 1: triangulation not valid; 2:
             // triangulation valid
      {
        auto polyCells = dataSetAsPD->GetPolys();
        auto lineCells = dataSetAsPD->GetLines();
        if(polyCells->GetNumberOfCells() > 0) {
          it = ttkAlgorithm::DataSetToTriangulationMap.find((void *)polyCells);
          if(it != ttkAlgorithm::DataSetToTriangulationMap.end())
            triangulationStatus
              = polyCells->GetMTime() != it->second.second ? 1 : 2;
        } else if(lineCells->GetNumberOfCells() > 0) {
          it = ttkAlgorithm::DataSetToTriangulationMap.find((void *)lineCells);
          if(it != ttkAlgorithm::DataSetToTriangulationMap.end())
            triangulationStatus
              = lineCells->GetMTime() != it->second.second ? 1 : 2;
        }
      }

      if(triangulationStatus == 2) {
        this->PrintMsg("Returning already initilized triangulation",
                       ttk::debug::Priority::DETAIL);
        return it->second.first;
      }
      if(triangulationStatus == 1) {
        this->PrintMsg("Chached triangulation no longer valid",
                       ttk::debug::Priority::DETAIL);
      }

      this->PrintMsg(
        "Initializing triangulation", 0, ttk::debug::Priority::DETAIL);
      auto newTriangulation = new ttk::Triangulation();

      // init points
      if(dataSet->GetNumberOfPoints() > 0) {
        auto points = dataSetAsPD->GetPoints();
        auto pointDataType = points->GetDataType();

        if(pointDataType != VTK_FLOAT && pointDataType != VTK_DOUBLE) {
          this->PrintErr("Unable to initialize 'ttk::Triangulation' for point "
                         "precision other than 'float' or 'double'.");
          delete newTriangulation;
          return 0;
        }

        newTriangulation->setInputPoints(dataSet->GetNumberOfPoints(),
                                         points->GetVoidPointer(0),
                                         pointDataType == VTK_DOUBLE);
      }

      // init cells
      {
        auto polyCells = dataSetAsPD->GetPolys();
        auto lineCells = dataSetAsPD->GetLines();
        if(polyCells->GetNumberOfCells() > 0) {
          // 2D
          newTriangulation->setInputCells(
            polyCells->GetNumberOfCells(), polyCells->GetPointer());

          ttkAlgorithm::DataSetToTriangulationMap.insert(
            {(void *)polyCells, {newTriangulation, polyCells->GetMTime()}});
        } else if(lineCells->GetNumberOfCells() > 0) {
          // 1D
          newTriangulation->setInputCells(
            lineCells->GetNumberOfCells(), lineCells->GetPointer());

          ttkAlgorithm::DataSetToTriangulationMap.insert(
            {(void *)lineCells, {newTriangulation, lineCells->GetMTime()}});
        }

        this->PrintMsg("Initializing Triangulation", 1,
                       ttk::debug::LineMode::REPLACE,
                       ttk::debug::Priority::DETAIL);
        return newTriangulation;
      }

      break;
    }

    // =====================================================================
    case VTK_IMAGE_DATA: {
      // auto dataSetAsID = (vtkImageData*) dataSet;

      // int extents[6];
      // dataSetAsID->GetExtent(extents);

      // double origin[3];
      // dataSetAsID->GetOrigin(origin);

      // double spacing[3];
      // dataSetAsID->GetSpacing(spacing);

      // int gridDimensions[3];
      // dataSetAsID->GetDimensions(gridDimensions);

      // double firstPoint[3];
      // firstPoint[0] = origin[0] + extents[0] * spacing[0];
      // firstPoint[1] = origin[1] + extents[2] * spacing[1];
      // firstPoint[2] = origin[2] + extents[4] * spacing[2];

      // newTriangulation->setInputGrid(
      //     firstPoint[0], firstPoint[1], firstPoint[2],
      //     spacing[0], spacing[1], spacing[2],
      //     gridDimensions[0], gridDimensions[1], gridDimensions[2]
      // );

      // cellKey = (void*)dataSet;
      // initTime = -1;
      break;
    }

    // =====================================================================
    // UNSUPPORTED DATA TYPE
    // =====================================================================
    default: {
    }
  }

  // if(cellKey==nullptr || initTime<0 || triangulation==nullptr){
  // return 0;
  // }

  this->PrintErr("Unable to get/create triangulation for '"
                 + std::string(dataSet->GetClassName()) + "'");

  return nullptr;
}

template <class vtkDataType>
int prepOutput(vtkInformation *info, std::string className) {
  auto output = vtkDataObject::GetData(info);
  if(!output || !output->IsA(className.data())) {
    auto newOutput = vtkSmartPointer<vtkDataType>::New();
    info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
  }
  return 1;
}

int ttkAlgorithm::RequestDataObject(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // for each output
  for(int i = 0; i < this->GetNumberOfOutputPorts(); ++i) {
    auto outInfo = outputVector->GetInformationObject(i);

    auto outputPortInfo = this->GetOutputPortInformation(i);

    // always request output type again for dynamic filter outputs
    this->FillOutputPortInformation(i, outputPortInfo);

    if(outputPortInfo->Has(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT())) {
      // Set output data type to input data type at specified port
      auto inPortIndex
        = outputPortInfo->Get(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());
      if(inPortIndex >= this->GetNumberOfInputPorts()) {
        this->PrintErr("Input port index " + std::to_string(inPortIndex)
                       + " specified by 'SAME_DATA_TYPE_AS_INPUT_PORT' key of "
                         "output port is out of range ("
                       + std::to_string(this->GetNumberOfInputPorts())
                       + " input ports).");
        return 0;
      }
      auto inInfo = inputVector[inPortIndex]->GetInformationObject(0);
      if(!inInfo) {
        this->PrintErr(
          "No information object at port " + std::to_string(inPortIndex)
          + " specified by 'SAME_DATA_TYPE_AS_INPUT_PORT' key of output port.");
        return 0;
      }

      auto input = vtkDataObject::GetData(inInfo);
      auto output = vtkDataObject::GetData(outInfo);

      if(!output || !output->IsA(input->GetClassName())) {
        auto newOutput
          = vtkSmartPointer<vtkDataObject>::Take(input->NewInstance());
        outputPortInfo->Set(
          vtkDataObject::DATA_TYPE_NAME(), input->GetClassName());
        outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
      }
    } else {
      // Explicitly create output by data type name
      if(!outputPortInfo->Has(vtkDataObject::DATA_TYPE_NAME())) {
        this->PrintErr("DATA_TYPE_NAME of output port " + std::to_string(i)
                       + " not specified");
        return 0;
      }
      std::string outputType
        = outputPortInfo->Get(vtkDataObject::DATA_TYPE_NAME());

      if(outputType == "vtkUnstructuredGrid") {
        prepOutput<vtkUnstructuredGrid>(outInfo, "vtkUnstructuredGrid");
      } else if(outputType == "vtkPolyData") {
        prepOutput<vtkPolyData>(outInfo, "vtkPolyData");
      } else if(outputType == "vtkMultiBlockDataSet") {
        prepOutput<vtkMultiBlockDataSet>(outInfo, "vtkMultiBlockDataSet");
      } else if(outputType == "vtkTable") {
        prepOutput<vtkTable>(outInfo, "vtkTable");
      } else if(outputType == "vtkImageData") {
        prepOutput<vtkImageData>(outInfo, "vtkImageData");
      } else {
        this->PrintErr("Unsupported data type for output[" + std::to_string(i)
                       + "]: " + outputType);
        return 0;
      }
    }

    this->PrintMsg(
      "Created '"
        + std::string(outputPortInfo->Get(vtkDataObject::DATA_TYPE_NAME()))
        + "' at output port " + std::to_string(i),
      ttk::debug::Priority::VERBOSE);
  }

  return 1;
}

//==============================================================================
int ttkAlgorithm::ProcessRequest(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  // 1. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_OBJECT())) {
    this->PrintMsg(
      "Processing REQUEST_DATA_OBJECT", ttk::debug::Priority::VERBOSE);
    return this->RequestDataObject(request, inputVector, outputVector);
  }

  // 2. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_INFORMATION())) {
    this->PrintMsg(
      "Processing REQUEST_INFORMATION", ttk::debug::Priority::VERBOSE);
    return this->RequestInformation(request, inputVector, outputVector);
  }

  // 3. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_TIME())) {
    this->PrintMsg(
      "Processing REQUEST_UPDATE_TIME", ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateTime(request, inputVector, outputVector);
  }

  // 4. Pass
  if(request->Has(
       vtkCompositeDataPipeline::REQUEST_TIME_DEPENDENT_INFORMATION())) {
    this->PrintMsg("Processing REQUEST_TIME_DEPENDENT_INFORMATION",
                   ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateTimeDependentInformation(
      request, inputVector, outputVector);
  }

  // 5. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_EXTENT())) {
    this->PrintMsg(
      "Processing REQUEST_UPDATE_EXTENT", ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  // 6. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_NOT_GENERATED())) {
    this->PrintMsg(
      "Processing REQUEST_DATA_NOT_GENERATED", ttk::debug::Priority::VERBOSE);
    return this->RequestDataNotGenerated(request, inputVector, outputVector);
  }

  // 7. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA())) {
    this->PrintMsg("Processing REQUEST_DATA", ttk::debug::Priority::VERBOSE);
    this->PrintMsg(ttk::debug::Separator::L0);
    return this->RequestData(request, inputVector, outputVector);
  }

  this->PrintErr("Unsupported pipeline pass:");
  request->Print(cout);

  return 0;
}