/// \ingroup vtk
/// \class ttkPathCompression
/// \author Robin Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief Compute the marching tetrahedra algorithm on a dataset
//
/// Input a dataset with an array containing labels, this application will
/// generate separating geometry between the labels.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkMarchingTetrahedra.h>

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
  std::string outputName{"MarchingTetOutput"};

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

  auto marchingTets = vtkSmartPointer<ttkMarchingTetrahedra>::New();

  vtkDataArray *defaultArray = nullptr;
  auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
  reader->SetFileName(inputFilePaths[0].data());
  reader->Update();

  auto inputDataObject = reader->GetOutput();
  if(!inputDataObject)
    return msg.printErr("Unable to read input file `" + inputFilePaths[0]
                        + "' :(");

  auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

  marchingTets->SetInputDataObject(0, reader->GetOutput());

  if(!defaultArray) {
    defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
    if(!defaultArray)
      defaultArray = inputAsVtkDataSet->GetCellData()->GetArray(0);
  }

  if(!inputArrayNames.size()) {
    if(defaultArray)
      inputArrayNames.emplace_back(defaultArray->GetName());
  }

  marchingTets->SetSurfaceMode(
    ttk::MarchingTetrahedra::SURFACE_MODE::SM_SEPARATORS);
  marchingTets->SetInputArrayToProcess(0, 0, 0, 0, inputArrayNames[0].data());
  marchingTets->Update();

  auto output = marchingTets->GetOutputDataObject(0);
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
