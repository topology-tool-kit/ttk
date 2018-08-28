#include                  <ttkManifoldLearning.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkManifoldLearning)

  int ttkManifoldLearning::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){
    Memory m;

    {
      stringstream msg;
      msg << "[ttkManifoldLearning] Memory usage: " << m.getElapsedUsage() 
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }
