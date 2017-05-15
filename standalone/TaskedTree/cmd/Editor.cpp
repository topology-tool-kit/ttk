/*
 * file:                  Editor.cpp
 * description:           Data-structures and processing.
 * author:                Your Name Here <Your Email Address Here>.
 * date:                  The Date Here.
 */

#include                  <Editor.h>

#include <vtksys/SystemTools.hxx>

Editor::Editor(){
   grid_       = NULL;
   lastObject_ = true;

   debug_       = 1;
   contourTree_ = vtkTaskedTree::New();

   core_           = -INT_MAX;
   fieldId_        = -INT_MAX;
   treeType_       = -INT_MAX;
   lessPartitions_ = false;
   partitionNum_   = -INT_MAX;
   threshold_      = -DBL_MAX;
   method_         = -INT_MAX;
}

Editor::~Editor(){

  // delete the mesh reader and the associated data at once
  contourTree_->Delete();
}

int Editor::execute(){

  contourTree_->setDebugLevel(debug_);
  contourTree_->SetdebugLevel_(debug_);
  contourTree_->SetThreadNumber(core_);
  contourTree_->SetUseAllCores(false);

  contourTree_->SetInputData(grid_);
  contourTree_->SetScalarFieldId(fieldId_);
  contourTree_->SettreeType_(treeType_);

  contourTree_->Update();

  grid_->ShallowCopy(contourTree_->GetOutput());

  return 0;
}

int Editor::init(int &argc, char **argv){

  CommandLineParser parser;

  // specify argument "-g" : The path of the grid
  parser.setArgument("g", &inputFilePath_, "Path to the input 3D grid");
  //parser.setIntArgument("d", &debug_, "Debug Level", true);
  parser.setArgument("n", &core_, "Thread number", true);
  parser.setArgument("f", &fieldId_, "Field identifier", true);
  parser.setArgument("s", &threshold_, "Simplification threshold between 0 and 1", true);
  parser.setArgument("m", &method_, "Simplification method : 0 persist, 1 Vertices ...", true);
  parser.setArgument("t", &treeType_, "type of tree : 2 is CT", true);

  // now parse the command line
  parser.parse(argc, argv);

  // set default values
  debug_ = ttk::globalDebugLevel_;

  if (debug_ == -INT_MAX) {
     debug_ = 1;
  }

  if (fieldId_ == -INT_MAX) {
     fieldId_ = 0;
  }

  if (core_ == -INT_MAX) {
     core_ = 1;
  }

  if (threshold_ == -DBL_MAX) {
     threshold_ = 0;
  }

  if (method_ == -INT_MAX) {
     method_ = 0;
  }

  if (treeType_ == -INT_MAX) {
     treeType_ = 2;
  }

  // now load the data to the editor
  loadData();

  return 0;
}

int Editor::loadData(){

  // create a reader object
  std::string extension =
          vtksys::SystemTools::GetFilenameLastExtension(inputFilePath_);

  if (extension == ".vtu"){
    grid_ = ReadAnXMLFile<vtkXMLUnstructuredGridReader> (inputFilePath_.c_str());
  } else if (extension == ".vti"){
    grid_ = ReadAnXMLFile<vtkXMLImageDataReader> (inputFilePath_.c_str());
  } else {
    cerr << "Bad format, need vtu" << endl;
    return -1;
  }

  {
    stringstream msg;
    msg << "[Editor]   done! (read "
      << grid_->GetNumberOfPoints()
      << " vertices, "
      << grid_->GetNumberOfCells()
      << " cells)" << endl;
    dMsg(cout, msg.str(), 1);
  }

  return 0;
}

int Editor::saveData(const string &fileName) const{

  vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();

  writer->SetFileName(fileName.data());
  writer->SetInputData(
          (vtkUnstructuredGrid*)contourTree_->GetOutput());
  writer->Write();

  writer->Delete();

  return 0;
}
