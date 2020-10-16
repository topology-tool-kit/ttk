/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for Reeb Space computation.

// include the local headers
#include <CommandLineParser.h>
#include <ttkReebSpace.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputPathPrefix{"output"};
  bool listArrays{false};
  bool forceOffset{false};

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

    parser.setOption("F", &forceOffset, "Force custom offset field (array #1)");

    parser.parse(argc, argv);
  }

  ttk::Debug msg;
  msg.setDebugMsgPrefix("ReebSpace");

  vtkNew<ttkReebSpace> rs{};

  vtkDataArray *defaultArray = nullptr;
  for(size_t i = 0; i < inputFilePaths.size(); i++) {
    // init a reader that can parse any vtkDataObject stored in xml format
    vtkNew<vtkXMLGenericDataObjectReader> reader{};
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
      // feed input object to the filter
      rs->SetInputDataObject(i, reader->GetOutput());

      // default arrays
      if(!defaultArray) {
        defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
        if(!defaultArray)
          defaultArray = inputAsVtkDataSet->GetCellData()->GetArray(0);
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
  if(!inputArrayNames.size()) {
    if(defaultArray)
      inputArrayNames.push_back(defaultArray->GetName());
  }
  for(size_t i = 0; i < inputArrayNames.size(); i++)
    rs->SetInputArrayToProcess(i, 0, 0, 0, inputArrayNames[i].data());

  // TODO manage the offset array?

  // ---------------------------------------------------------------------------
  // Execute the filter
  // ---------------------------------------------------------------------------
  rs->SetForceInputOffsetScalarField(forceOffset);
  rs->Update();

  // ---------------------------------------------------------------------------
  // If output prefix is specified then write all output objects to disk
  // ---------------------------------------------------------------------------
  if(!outputPathPrefix.empty()) {
    for(int i = 0; i < rs->GetNumberOfOutputPorts(); i++) {
      auto output = rs->GetOutputDataObject(i);
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
