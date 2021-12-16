/// TODO 12: Add your information
/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkCompactTriangulationPreconditioning.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program variables
  // ---------------------------------------------------------------------------
  int bucketThreshold = 500;
  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputPathPrefix{"output"};
  bool listArrays{false};

  // ---------------------------------------------------------------------------
  // Set program variables based on command line arguments
  // ---------------------------------------------------------------------------
  {
    ttk::CommandLineParser parser;

    // -------------------------------------------------------------------------
    // Standard options and arguments
    // -------------------------------------------------------------------------
    parser.setArgument(
      "i", &inputFilePaths, "Input data-sets (*.vti, *vtu, *vtp)", false);
    parser.setArgument("a", &inputArrayNames, "Input array names", true);
    parser.setArgument(
      "o", &outputPathPrefix, "Output file prefix (no extension)", true);
    parser.setOption("l", &listArrays, "List available arrays");

    // -------------------------------------------------------------------------
    // TODO 13: Declare custom arguments and options
    // -------------------------------------------------------------------------
    parser.setArgument("b", &bucketThreshold, "Bucket threshold", true);

    parser.parse(argc, argv);
  }

  // ---------------------------------------------------------------------------
  // Command line output messages.
  // ---------------------------------------------------------------------------
  ttk::Debug msg;
  msg.setDebugMsgPrefix("CompactTriangulationPreconditioning");

  // ---------------------------------------------------------------------------
  // Initialize ttkCompactTriangulationPreconditioning module (adjust
  // parameters)
  // ---------------------------------------------------------------------------
  auto compactTriangulationPreconditioning
    = vtkSmartPointer<ttkCompactTriangulationPreconditioning>::New();

  // ---------------------------------------------------------------------------
  // TODO 14: Pass custom arguments and options to the module
  // ---------------------------------------------------------------------------
  // compactTriangulationPreconditioning->SetOutputArrayName(outputArrayName);
  compactTriangulationPreconditioning->SetThreshold(bucketThreshold);

  // ---------------------------------------------------------------------------
  // Read input vtkDataObjects (optionally: print available arrays)
  // ---------------------------------------------------------------------------
  for(size_t i = 0; i < inputFilePaths.size(); i++) {
    // init a reader that can parse any vtkDataObject stored in xml format
    auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(inputFilePaths[i].data());
    reader->Update();

    // check if input vtkDataObject was successfully read
    auto inputDataObject = reader->GetOutput();
    if(!inputDataObject) {
      msg.printErr("Unable to read input file `" + inputFilePaths[i] + "' :(");
      return 1;
    }

    auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

    // if requested print list of arrays, otherwise proceed with execution
    if(listArrays) {
      msg.printMsg(inputFilePaths[i] + ":");
      if(inputAsVtkDataSet) {
        // Point Data
        msg.printMsg("  PointData:");
        auto pointData = inputAsVtkDataSet->GetPointData();
        for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(pointData->GetArrayName(j)));

        // Cell Data
        msg.printMsg("  CellData:");
        auto cellData = inputAsVtkDataSet->GetCellData();
        for(int j = 0; j < cellData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(cellData->GetArrayName(j)));
      } else {
        msg.printErr("Unable to list arrays on file `" + inputFilePaths[i]
                     + "'");
        return 1;
      }
    } else {
      // feed input object to ttkCompactTriangulationPreconditioning filter
      compactTriangulationPreconditioning->SetInputDataObject(
        i, reader->GetOutput());
      auto pointData = inputAsVtkDataSet->GetPointData();
      for(int j = 0; j < pointData->GetNumberOfArrays(); j++) {
        inputArrayNames.push_back(pointData->GetArrayName(j));
      }
      auto cellData = inputAsVtkDataSet->GetCellData();
      for(int j = 0; j < cellData->GetNumberOfArrays(); j++) {
        inputArrayNames.push_back(cellData->GetArrayName(j));
      }
    }
  }

  // terminate program if it was just asked to list arrays
  if(listArrays) {
    return 0;
  }

  // ---------------------------------------------------------------------------
  // Specify which arrays of the input vtkDataObjects will be processed
  // ---------------------------------------------------------------------------
  vtkNew<vtkDataArraySelection> arraySelection;
  for(size_t i = 0; i < inputArrayNames.size(); i++) {
    arraySelection->EnableArray(inputArrayNames[i].data());
    compactTriangulationPreconditioning->SetDataArraySelection(arraySelection);
  }

  // ---------------------------------------------------------------------------
  // Execute ttkCompactTriangulationPreconditioning filter
  // ---------------------------------------------------------------------------
  compactTriangulationPreconditioning->Update();

  // ---------------------------------------------------------------------------
  // If output prefix is specified then write all output objects to disk
  // ---------------------------------------------------------------------------
  if(!outputPathPrefix.empty()) {
    for(int i = 0;
        i < compactTriangulationPreconditioning->GetNumberOfOutputPorts();
        i++) {
      auto output = compactTriangulationPreconditioning->GetOutputDataObject(i);
      auto writer
        = vtkXMLDataObjectWriter::NewWriter(output->GetDataObjectType());

      std::string outputFileName = outputPathPrefix + "_port_"
                                   + std::to_string(i) + "."
                                   + writer->GetDefaultFileExtension();
      msg.printMsg("Writing output file `" + outputFileName + "'...");
      writer->SetInputDataObject(output);
      writer->SetFileName(outputFileName.data());
      writer->Update();
    }
  }

  return 0;
}
