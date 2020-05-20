/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

#include <vtkSmartPointer.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtksys/CommandLineArguments.hxx>

#include <ttkHelloWorld.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program Variables
  // ---------------------------------------------------------------------------
  std::string inputFileName{""};
  std::string inputScalarFieldName{""};
  int debugLevel{0};
  int threadNumber{1};
  std::string outputFileName;

  // ---------------------------------------------------------------------------
  // Set Program Variables based on Command Line Arguments
  // ---------------------------------------------------------------------------
  {
    // init command line parser
    vtksys::CommandLineArguments arg;
    arg.StoreUnusedArguments(1);
    arg.Initialize(argc, argv);

    // specify arguments
    typedef vtksys::CommandLineArguments argT;
    arg.AddArgument(
      "-i", argT::SPACE_ARGUMENT, &inputFileName, "Path to input vtkDataSet");
    arg.AddArgument(
      "-o", argT::SPACE_ARGUMENT, &outputFileName, "Path to output vtkDataSet");
    arg.AddArgument(
      "-f", argT::SPACE_ARGUMENT, &inputScalarFieldName, "Scalar field name");
    arg.AddArgument(
      "-d", argT::SPACE_ARGUMENT, &debugLevel, "Debug Level");
    arg.AddArgument(
      "-t", argT::SPACE_ARGUMENT, &threadNumber, "Thread Number");

    // check if required arguments are set correctly
    if(!arg.Parse() || inputFileName.empty() || inputScalarFieldName.empty()
       || outputFileName.empty()) {
      std::cout << "Error parsing arguments" << std::endl;
      std::cout << arg.GetHelp() << std::endl;
      return 0;
    }
  }

  // ---------------------------------------------------------------------------
  // Execute Test
  // ---------------------------------------------------------------------------
  {
    ttk::globalDebugLevel_ = debugLevel;

    auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(inputFileName.data());

    auto helloWorld = vtkSmartPointer<ttkHelloWorld>::New();
    helloWorld->SetInputConnection(0, reader->GetOutputPort(0));
    helloWorld->SetInputArrayToProcess(0, 0, 0, 0, inputScalarFieldName.data());
    helloWorld->SetUseAllCores(false);
    helloWorld->SetThreadNumber(threadNumber);

    auto writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
    writer->SetInputConnection(0, helloWorld->GetOutputPort(0));
    writer->SetFileName(outputFileName.data());
    writer->Update();
  }

  return 1;
}
