/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>.
/// \date February 2017.
///
/// \brief Command line program for critical point computation.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkMergeAndContourTree.h>
#include <ttkMergeTree.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program variables
  // ---------------------------------------------------------------------------
  std::string inputFilePath = "";
  std::string inputArrayName = "";
  bool listArrays{false};
  bool compare{false};
  int threadNumber{112};
  int repetitions{1};

  // ---------------------------------------------------------------------------
  // Set program variables based on command line arguments
  // ---------------------------------------------------------------------------
  {
    ttk::CommandLineParser parser;

    // -------------------------------------------------------------------------
    // Standard options and arguments
    // -------------------------------------------------------------------------
    parser.setArgument(
      "i", &inputFilePath, "Input data-set (*.vti, *vtu, *vtp)", false);
    parser.setArgument("a", &inputArrayName, "Input array name", false);
    parser.setArgument("r", &repetitions,
                       "Number of times you want to run the algorithms", true);
    parser.setOption("l", &listArrays, "List available arrays");
    parser.setOption("c", &compare, "Also build the TTK FTM Tree to compare");
    parser.setArgument(
      "n", &threadNumber, "The number of OMP Threads to run the filters", true);

    parser.parse(argc, argv);
  }

  // ---------------------------------------------------------------------------
  // Command line output messages.
  // ---------------------------------------------------------------------------
  ttk::Debug msg;
  msg.setDebugMsgPrefix("MergeTree");

  // ---------------------------------------------------------------------------
  // Initialize ttkPairExtrema module (adjust parameters)
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Read input vtkDataObjects (optionally: print available arrays)
  // ---------------------------------------------------------------------------
  vtkDataArray *defaultArray = nullptr;
  // init a reader that can parse any vtkDataObject stored in xml format
  auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
  reader->SetFileName(inputFilePath.data());
  reader->Update();

  // check if input vtkDataObject was successfully read
  auto inputDataObject = reader->GetOutput();
  if(!inputDataObject) {
    msg.printErr("Unable to read input file `" + inputFilePath + "' :(");
    return 1;
  }

  auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

  // if requested print list of arrays, otherwise proceed with execution
  if(listArrays) {
    msg.printMsg(inputFilePath + ":");
    if(inputAsVtkDataSet) {
      // Point Data
      msg.printMsg("  PointData:");
      auto pointData = inputAsVtkDataSet->GetPointData();
      for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
        msg.printMsg("    - " + std::string(pointData->GetArrayName(j)));
      return 0;
    } else {
      msg.printErr("Unable to list arrays on file `" + inputFilePath + "'");
      return 1;
    }
  }

  // ---------------------------------------------------------------------------
  // Specify which arrays of the input vtkDataObjects will be processed
  // ---------------------------------------------------------------------------
  if(!defaultArray) {
    defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
  }
  if(inputArrayName == "") {
    if(defaultArray)
      inputArrayName = defaultArray->GetName();
  }

  for(int i = 0; i < repetitions; i++) {
    auto mergeTree = vtkSmartPointer<ttkMergeTree>::New();

    mergeTree->SetInputDataObject(0, inputDataObject);

    mergeTree->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());

    mergeTree->SetUseAllCores(false);
    mergeTree->SetThreadNumber(threadNumber);
    mergeTree->Modified();
    mergeTree->Update();
  }

  if(compare) {
    for(int i = 0; i < repetitions; i++) {
      msg.setDebugMsgPrefix("FTMTree");
      auto contourTree = vtkSmartPointer<ttkMergeAndContourTree>::New();
      contourTree->SetInputDataObject(0, reader->GetOutput());
      contourTree->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());
      contourTree->SetUseAllCores(false);
      contourTree->SetThreadNumber(threadNumber);
      contourTree->SetTreeType(2); // Contour tree
      contourTree->Modified();
      contourTree->Update();
    }
  }

  return 0;
}
