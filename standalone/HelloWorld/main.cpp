/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkHelloWorld.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program variables
  // ---------------------------------------------------------------------------
  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputPathPrefix{""};
  int listArrays{0};

  // ---------------------------------------------------------------------------
  // Set program variables based on command line arguments
  // ---------------------------------------------------------------------------
  {
    ttk::CommandLineParser parser;
    parser.setArgument(
      "i", &inputFilePaths, "Input data-sets (*.vti, *vtu, *vtp)", false);
    parser.setArgument("a", &inputArrayNames, "Input array names", true);
    parser.setArgument(
      "o", &outputPathPrefix, "Output file prefix (no extension)", true);
    parser.setArgument("l", &listArrays, "List available arrays", true);

    parser.parse(argc, argv);
    parser.printArgs();
  }

  // ---------------------------------------------------------------------------
  // Initialize ttkHelloWorld module (adjust parameters)
  // ---------------------------------------------------------------------------
  auto helloWorld = vtkSmartPointer<ttkHelloWorld>::New();
  // helloWorld->SetOutputArrayName("AveragedArray");

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
      std::cout << "ERROR: Unable to read input vtkDataObject " << i << "\n";
      return 0;
    }

    // if requested print list of arrays, otherwise proceed with execution
    if(listArrays) {
      std::cout << "vtkDataObject " << i << ":\n";
      if(auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject)) {
        // Point Data
        std::cout << "  PointData:\n";
        auto pointData = inputAsVtkDataSet->GetPointData();
        for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
          std::cout << "    " << pointData->GetArrayName(j) << "\n";

        // Cell Data
        std::cout << "  CellData:\n";
        auto cellData = inputAsVtkDataSet->GetCellData();
        for(int j = 0; j < cellData->GetNumberOfArrays(); j++)
          std::cout << "    " << cellData->GetArrayName(j) << "\n";
      } else {
        std::cout << "ERROR: Unable to list arrays of non vtkDataSet input.";
        return 0;
      }
    } else {
      // feed input object to ttkHelloWorld filter
      helloWorld->SetInputDataObject(i, reader->GetOutput());
    }
  }

  // terminate program if it was just asked to list arrays
  if(listArrays) {
    return 1;
  }

  // ---------------------------------------------------------------------------
  // Specify which arrays of the input vtkDataObjects will be processed
  // ---------------------------------------------------------------------------
  for(size_t i = 0; i < inputArrayNames.size(); i++)
    helloWorld->SetInputArrayToProcess(i, 0, 0, 0, inputArrayNames[i].data());

  // ---------------------------------------------------------------------------
  // Execute ttkHelloWorld filter
  // ---------------------------------------------------------------------------
  helloWorld->Update();

  // ---------------------------------------------------------------------------
  // If output prefix is specified then write all output objects to disk
  // ---------------------------------------------------------------------------
  if(!outputPathPrefix.empty()) {
    for(int i = 0; i < helloWorld->GetNumberOfOutputPorts(); i++) {
      auto output = helloWorld->GetOutputDataObject(i);
      auto writer
        = vtkXMLDataObjectWriter::NewWriter(output->GetDataObjectType());
      writer->SetInputDataObject(output);
      writer->SetFileName((outputPathPrefix + "_port_" + std::to_string(i) + "."
                           + writer->GetDefaultFileExtension())
                            .data());
      writer->Update();
    }
  }

  return 1;
}