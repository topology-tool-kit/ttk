/// \ingroup vtk
/// \class ttkPathCompression
/// \author Robin Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief Compute the ascending, descending, and MS segmentation hash from a
/// dataset.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkPathCompression.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputName{"MSSegmentation"};

  {
    ttk::CommandLineParser parser;

    parser.setArgument(
      "i", &inputFilePaths, "Input data-sets (*.vti, *vtu, *vtp)", false);
    parser.setArgument("a", &inputArrayNames, "Input array name", true);

    parser.setArgument(
      "o", &outputName, "Output file name (no extension)", true);

    parser.parse(argc, argv);
  }

  ttk::Debug msg;
  msg.setDebugMsgPrefix("PathCompression");

  auto pathCompression = vtkSmartPointer<ttkPathCompression>::New();

  vtkDataArray *defaultArray = nullptr;
  auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
  reader->SetFileName(inputFilePaths[0].data());
  reader->Update();

  auto inputDataObject = reader->GetOutput();
  if(!inputDataObject) {
    msg.printErr("Unable to read input file `" + inputFilePaths[0] + "' :(");
    return 1;
  }

  auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

  pathCompression->SetInputDataObject(0, reader->GetOutput());

  if(!defaultArray) {
    defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
    if(!defaultArray)
      defaultArray = inputAsVtkDataSet->GetCellData()->GetArray(0);
  }

  if(!inputArrayNames.size()) {
    if(defaultArray)
      inputArrayNames.emplace_back(defaultArray->GetName());
  }

  pathCompression->SetComputeAscendingSegmentation(true);
  pathCompression->SetComputeDescendingSegmentation(true);
  pathCompression->SetComputeMSSegmentationHash(true);
  pathCompression->SetInputArrayToProcess(
    0, 0, 0, 0, inputArrayNames[0].data());
  pathCompression->Update();

  auto output = pathCompression->GetOutputDataObject(0);
  auto writer = vtkSmartPointer<vtkXMLWriter>::Take(
    vtkXMLDataObjectWriter::NewWriter(output->GetDataObjectType()));

  std::string outputFileName
    = outputName + "." + writer->GetDefaultFileExtension();
  msg.printMsg("Writing output file `" + outputFileName + "'...");
  writer->SetInputDataObject(output);
  writer->SetFileName(outputFileName.data());
  writer->Update();

  return 0;
}
