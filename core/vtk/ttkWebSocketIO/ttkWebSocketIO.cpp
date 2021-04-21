#include <ttkWebSocketIO.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>

#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkWebSocketIO);

ttkWebSocketIO::ttkWebSocketIO() {
  this->LastInput = vtkSmartPointer<vtkUnstructuredGrid>::New();
  this->LastOutput = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkWebSocketIO::~ttkWebSocketIO() {
}

int ttkWebSocketIO::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkWebSocketIO::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkWebSocketIO::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  const bool neededUpdate = this->GetNeedsUpdate();
  this->SetNeedsUpdate(false);

  // close existing server if it has a different port
  if(this->isListening() && (this->getPortNumber() != this->PortNumber))
    this->stopServer();

  try {
    if(!this->isListening())
      this->startServer(this->PortNumber);
  } catch(const std::exception &e) {
    this->printErr("Unable to start server:\n" + std::string(e.what()));
    return 0;
  }

  if(neededUpdate) {
    // if needed update then the client send data to the server which was parsed
    // and is now passed on as the output of the filter
    auto output = vtkMultiBlockDataSet::GetData(outputVector);
    output->ShallowCopy(this->LastOutput);
  } else {
    // in any other case (i.e., if the input or parameters changed) then send
    // the last input to the client
    auto input = vtkDataObject::GetData(inputVector[0]);
    this->LastInput
      = vtkSmartPointer<vtkDataObject>::Take(input->NewInstance());
    this->LastInput->ShallowCopy(input);
    if(!this->processEvent("on_message", "RequestInputVtkDataSet"))
      return 0;
  }
  return 1;
}

int ttkWebSocketIO::processEvent(const std::string &eventName,
                                 const std::string &eventData) {
  if(eventData.compare("RequestInputVtkDataSet") == 0) {
    if(!this->SendVtkDataObject(this->LastInput))
      return 0;
  } else if(eventData.rfind("{\"vtkDataSet", 12) == 0) {
    if(!this->ParseVtkDataObjectFromJSON(eventData))
      return 0;
  }

  return WebSocketIO::processEvent(eventName, eventData);
}

template <typename T>
T jsonGetValue(const boost::property_tree::ptree &pt,
               const boost::property_tree::ptree::key_type &key) {
  return pt.get_child(key).get_value<T>();
}

template <typename T>
void jsonArrayToArray(const boost::property_tree::ptree &pt,
                      const boost::property_tree::ptree::key_type &key,
                      T *result) {
  std::vector<double> values;
  ttkUtils::stringListToDoubleVector(
    jsonGetValue<std::string>(pt, key), values);

  for(size_t i = 0, j = values.size(); i < j; i++) {
    result[i] = (T)values[i];
  }
}

bool jsonHasChild(const boost::property_tree::ptree &pt,
                  const boost::property_tree::ptree::key_type &key) {
  return pt.find(key) != pt.not_found();
}

int ttkWebSocketIO::ParseVtkDataObjectFromJSON(const std::string &json) {
  ttk::Timer timer;
  this->printMsg("Parsing vtkDataObject from JSON string", 0, 0,
                 ttk::debug::LineMode::REPLACE);

  this->LastOutput = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  // parse lastClientInput into json
  boost::property_tree::ptree pt;
  try {
    std::stringstream ss;
    ss << json;
    boost::property_tree::read_json(ss, pt);
  } catch(const std::exception &e) {
    this->printErr("Unable to parse JSON message: " + json);
    return 0;
  }

  // check if a path exists in property_tree or a key exists in json object

  size_t b = 0;
  for(auto it = pt.begin(); it != pt.end(); ++it) {
    auto block = vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto jsonVtkDataSet = it->second;
    if(!jsonHasChild(jsonVtkDataSet, "className")
       || !jsonHasChild(jsonVtkDataSet, "points")
       || !jsonHasChild(jsonVtkDataSet, "cells")
       || !jsonHasChild(jsonVtkDataSet, "pointData")
       || !jsonHasChild(jsonVtkDataSet, "cellData")
       || !jsonHasChild(jsonVtkDataSet, "fieldData")) {
      this->printErr("Invalid vtkDataSet serialization.");
      return 0;
    }

    if(jsonGetValue<std::string>(jsonVtkDataSet, "className")
         .compare("vtkUnstructuredGrid")
       != 0) {
      this->printErr(
        "Currently WebSocketIO only supports parsing 'vtkUnstructuredGrids'.");
      return 0;
    }

    // parse points
    {
      auto nPoints
        = jsonGetValue<int>(jsonVtkDataSet, "points.coordinates.nTuples");

      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(nPoints);

      jsonArrayToArray<float>(jsonVtkDataSet, "points.coordinates.data",
                              (float *)ttkUtils::GetVoidPointer(points));

      block->SetPoints(points);
    }

    // parse cells
    {
      // VTK CELL TYPES MAP
      constexpr int vtkCellsTypeHash[20] = {VTK_EMPTY_CELL,
                                            VTK_VERTEX,
                                            VTK_LINE,
                                            VTK_TRIANGLE,
                                            VTK_TETRA,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_VOXEL,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET,
                                            VTK_CONVEX_POINT_SET};

      auto nOffsets
        = jsonGetValue<int>(jsonVtkDataSet, "cells.offsetsArray.nTuples");
      auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
      offsets->SetNumberOfTuples(nOffsets);
      auto offsetsData = (vtkIdType *)ttkUtils::GetVoidPointer(offsets);
      jsonArrayToArray<vtkIdType>(
        jsonVtkDataSet, "cells.offsetsArray.data", offsetsData);

      auto nConnectivity
        = jsonGetValue<int>(jsonVtkDataSet, "cells.connectivityArray.nTuples");
      auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
      connectivity->SetNumberOfTuples(nConnectivity);
      jsonArrayToArray<vtkIdType>(
        jsonVtkDataSet, "cells.connectivityArray.data",
        (vtkIdType *)ttkUtils::GetVoidPointer(connectivity));

      auto cells = vtkSmartPointer<vtkCellArray>::New();
      cells->SetData(offsets, connectivity);

      auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();
      cellTypes->SetNumberOfTuples(nOffsets - 1);
      auto cellTypesData = (unsigned char *)ttkUtils::GetVoidPointer(cellTypes);
      for(int i = 1; i < nOffsets; i++) {
        cellTypesData[i - 1]
          = vtkCellsTypeHash[offsetsData[i] - offsetsData[i - 1]];
      }

      block->SetCells(cellTypes, cells);
    }

    auto addArraysToAssociation = [](vtkFieldData *fd,
                                     const boost::property_tree::ptree &fdJSON,
                                     const std::string &fdKey) {
      // parse point data
      for(const auto &pda : fdJSON.get_child(fdKey)) {
        const auto &jsonArray = pda.second;
        size_t nComponents = jsonGetValue<size_t>(jsonArray, "nComponents");
        size_t nTuples = jsonGetValue<size_t>(jsonArray, "nTuples");

        auto array = vtkSmartPointer<vtkAbstractArray>::Take(
          vtkDataArray::CreateArray(jsonGetValue<int>(jsonArray, "dataType")));

        array->SetName(jsonGetValue<std::string>(jsonArray, "name").data());
        array->SetNumberOfComponents(nComponents);
        array->SetNumberOfTuples(nTuples);

        if(auto dataArray = vtkDataArray::SafeDownCast(array)) {
          switch(dataArray->GetDataType()) {
            vtkTemplateMacro(jsonArrayToArray<VTK_TT>(
              jsonArray, "data",
              (VTK_TT *)ttkUtils::GetVoidPointer(dataArray)));
          }
        } else if(auto stringArray = vtkStringArray::SafeDownCast(array)) {
          std::vector<std::string> values;
          ttkUtils::stringListToVector(
            jsonGetValue<std::string>(jsonArray, "data"), values);

          for(size_t i = 0; i < values.size(); i++)
            stringArray->SetValue(i, values[i]);
        } else {
          return 0;
        }

        fd->AddArray(array);
      }
      return 1;
    };

    if(!addArraysToAssociation(
         block->GetPointData(), jsonVtkDataSet, "pointData")) {
      this->printErr("Unsupported data type.");
      return 0;
    }
    if(!addArraysToAssociation(
         block->GetCellData(), jsonVtkDataSet, "cellData")) {
      this->printErr("Unsupported data type.");
      return 0;
    }
    if(!addArraysToAssociation(
         block->GetFieldData(), jsonVtkDataSet, "fieldData")) {
      this->printErr("Unsupported data type.");
      return 0;
    }

    this->LastOutput->SetBlock(b++, block);
  }

  this->printMsg(
    "Parsing vtkDataObject from JSON string", 1, timer.getElapsedTime());

  this->SetNeedsUpdate(true);

  return 1;
}

int ttkWebSocketIO::SendVtkDataObject(vtkDataObject *object) {
  ttk::Timer timer;
  this->printMsg(
    "Serializing vtkDataObject", 0, 0, ttk::debug::LineMode::REPLACE);

  // clear queue
  this->clearMessageQueue();

  auto objectAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(object->IsA("vtkMultiBlockDataSet"))
    objectAsMB->ShallowCopy(object);
  else
    objectAsMB->SetBlock(0, object);

  const size_t nBlocks = objectAsMB->GetNumberOfBlocks();
  for(size_t b = 0; b < nBlocks; b++) {
    auto block = vtkDataSet::SafeDownCast(objectAsMB->GetBlock(b));
    if(!block) {
      this->printErr("Input `vtkDataObject` is not a `vtkDataSet` or a "
                     "`vtkMultiBlockDataSet` containing only `vtkDataSets`.");
      return 0;
    }

    // send class name
    this->queueMessage("{ \"className\": \""
                       + std::string(block->GetClassName()) + "\" }");

    this->printMsg(
      "Serializing vtkDataObject", 0.2, 0, ttk::debug::LineMode::REPLACE);

    // send image meta data if last input is a vtkImageData
    if(block->IsA("vtkImageData")) {
      auto lastInputAsID = vtkImageData::SafeDownCast(block);
      int dimension[3];
      double origin[3];
      double spacing[3];
      lastInputAsID->GetDimensions(dimension);
      lastInputAsID->GetOrigin(origin);
      lastInputAsID->GetSpacing(spacing);

      this->queueMessage(
        "{"
        "\"dimension\": ["
        + std::to_string(dimension[0]) + "," + std::to_string(dimension[1])
        + "," + std::to_string(dimension[2])
        + "],"
          "\"origin\": ["
        + std::to_string(origin[0]) + "," + std::to_string(origin[1]) + ","
        + std::to_string(origin[2])
        + "],"
          "\"spacing\": ["
        + std::to_string(spacing[0]) + "," + std::to_string(spacing[1]) + ","
        + std::to_string(spacing[2])
        + "]"
          "}");
    }

    this->printMsg(
      "Serializing vtkDataObject", 0.4, 0, ttk::debug::LineMode::REPLACE);

    // send points if last input is a vtkPointSet
    auto blockAsPS = vtkPointSet::SafeDownCast(block);
    if(blockAsPS != nullptr) {
      auto points = blockAsPS->GetPoints();
      this->queueMessage(
        "{"
        "\"target\": \"points\","
        "\"name\": \"coordinates\","
        "\"dataType\":"
        + std::to_string(points ? points->GetDataType() : VTK_FLOAT)
        + ","
          "\"nTuples\":"
        + std::to_string(blockAsPS->GetNumberOfPoints())
        + ","
          "\"nComponents\": 3"
          "}");
      if(blockAsPS->GetNumberOfPoints() > 0 && points != nullptr)
        switch(points->GetDataType()) {
          vtkTemplateMacro(this->queueMessage(
            blockAsPS->GetNumberOfPoints() * sizeof(VTK_TT) * 3,
            ttkUtils::GetVoidPointer(points)));
        }
    }

    this->printMsg(
      "Serializing vtkDataObject", 0.6, 0, ttk::debug::LineMode::REPLACE);

    // send cells if last input is vtkUnstructuredGrid or vtkPolyData
    {
      vtkCellArray *cells = nullptr;

      if(block->IsA("vtkUnstructuredGrid")) {
        auto blockAsUG = vtkUnstructuredGrid::SafeDownCast(block);
        cells = blockAsUG->GetCells();
      } else if(block->IsA("vtkPolyData")) {
        auto blockAsPD = vtkPolyData::SafeDownCast(block);
        cells = blockAsPD->GetPolys() ? blockAsPD->GetPolys()
                                      : blockAsPD->GetLines();
      }

      if(cells) {
        {
          auto connectivityArray = cells->GetConnectivityArray();
          this->queueMessage(
            "{"
            "\"target\": \"cells\","
            "\"name\": \"connectivityArray\","
            "\"dataType\":"
            + std::to_string(VTK_ID_TYPE)
            + ","
              "\"nTuples\":"
            + std::to_string(connectivityArray->GetNumberOfTuples())
            + ","
              "\"nComponents\": 1"
              "}");
          if(connectivityArray->GetNumberOfTuples() > 0)
            this->queueMessage(
              connectivityArray->GetNumberOfTuples() * sizeof(vtkIdType),
              ttkUtils::GetVoidPointer(connectivityArray));
        }

        {
          auto offsetsArray = cells->GetOffsetsArray();
          this->queueMessage("{"
                             "\"target\": \"cells\","
                             "\"name\": \"offsetsArray\","
                             "\"dataType\":"
                             + std::to_string(VTK_ID_TYPE)
                             + ","
                               "\"nTuples\":"
                             + std::to_string(offsetsArray->GetNumberOfTuples())
                             + ","
                               "\"nComponents\": 1"
                               "}");
          if(offsetsArray->GetNumberOfTuples() > 0)
            this->queueMessage(
              offsetsArray->GetNumberOfTuples() * sizeof(vtkIdType),
              ttkUtils::GetVoidPointer(offsetsArray));
        }
      }
    }

    this->printMsg(
      "Serializing vtkDataObject", 0.8, 0, ttk::debug::LineMode::REPLACE);

    // send point, cell, and field data
    {
      std::vector<std::pair<std::string, vtkFieldData *>> attributes
        = {{"pointData", block->GetPointData()},
           {"cellData", block->GetCellData()},
           {"fieldData", block->GetFieldData()}};

      for(auto &attribute : attributes) {
        size_t nArrays = attribute.second->GetNumberOfArrays();

        for(size_t i = 0; i < nArrays; i++) {
          auto array = attribute.second->GetAbstractArray(i);

          if(array->IsA("vtkDataArray")) {
            auto dataArray = vtkDataArray::SafeDownCast(array);
            this->queueMessage(
              "{"
              "\"target\":\""
              + attribute.first
              + "\","
                "\"name\":\""
              + std::string(dataArray->GetName())
              + "\","
                "\"dataType\":"
              + std::to_string(dataArray->GetDataType())
              + ","
                "\"nTuples\":"
              + std::to_string(dataArray->GetNumberOfTuples())
              + ","
                "\"nComponents\":"
              + std::to_string(dataArray->GetNumberOfComponents()) + "}");
            if(dataArray->GetNumberOfValues() > 0)
              switch(dataArray->GetDataType()) {
                vtkTemplateMacro(this->queueMessage(
                  dataArray->GetNumberOfValues() * sizeof(VTK_TT),
                  ttkUtils::GetVoidPointer(dataArray)));
              }
          } else if(array->IsA("vtkStringArray")) {
            auto stringArray = vtkStringArray::SafeDownCast(array);
            std::string values;
            values += stringArray->GetValue(0);
            for(int j = 1; j < stringArray->GetNumberOfValues(); j++)
              values += "," + stringArray->GetValue(j);

            this->queueMessage(
              "{"
              "\"target\":\""
              + attribute.first
              + "\","
                "\"name\":\""
              + std::string(stringArray->GetName())
              + "\","
                "\"dataType\":"
              + std::to_string(stringArray->GetDataType())
              + ","
                "\"nTuples\":"
              + std::to_string(stringArray->GetNumberOfTuples())
              + ","
                "\"nComponents\":"
              + std::to_string(stringArray->GetNumberOfComponents())
              + ","
                "\"data\":\""
              + boost::replace_all_copy(values, "\"", "\\\"")
              + "\""
                "}");
          }
        }
      }
    }
  }

  this->printMsg("Serializing vtkDataObject", 1, timer.getElapsedTime());

  this->processMessageQueue();

  return 1;
}
