/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include <ttkLDistance.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkLDistance> program;

constexpr unsigned int str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

class myKeyHandler : public ttkKeyHandler {

public:
  int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key) {

    stringstream msg;
    msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
    dMsg(cout, msg.str(), infoMsg);

    // depending on the value of "key", trigger the right functions on the
    // program object (or its contained ttkObject_).
    switch(str2int(key.c_str())) {
      case str2int("space"):
        program.run();
        break;

      default:
        break;
    }

    return 0;
  }
};

int main(int argc, char **argv) {

  string distanceType;
  string distanceName = "";
  string dataset1 = "";
  string dataset2 = "";
  int id1 = -1;
  int id2 = -1;

  program.parser_.setArgument(
    "distance", &distanceType, "Distance type (1, 2, ..., inf)");
  program.parser_.setArgument("n1", &dataset1, "u field name", true);
  program.parser_.setArgument("n2", &dataset2, "v field name", true);
  program.parser_.setArgument(
    "fieldName", &distanceName, "output distance field name", true);

  program.parser_.setArgument("id1", &id1, "field id 1", true);
  program.parser_.setArgument("id2", &id2, "field id 2", true);

  int ret = program.init(argc, argv);
  if(ret != 0)
    return ret;

  // Change here arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetDistanceType(distanceType);

  // Check first scalar field.
  if(dataset1 != "")
    program.ttkObject_->SetScalarField1(dataset1);
  else if(id1 > -1)
    program.ttkObject_->SetScalarFieldId1(id1);
  else {
    stringstream msg;
    msg << "[MainThread] Badly specified scalar field 1." << endl;
    cerr << msg.str();
    return -1;
  }

  // Check second scalar field.
  if(dataset2 != "")
    program.ttkObject_->SetScalarField2(dataset2);
  else if(id2 > -1)
    program.ttkObject_->SetScalarFieldId2(id2);
  else {
    stringstream msg;
    msg << "[MainThread] Badly specified scalar field 2." << endl;
    cerr << msg.str();
    return -1;
  }

  if(distanceName != "")
    program.ttkObject_->SetDistanceFieldName(distanceName);

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler);

  // specify if the display of certain outputs should be disabled
  // vector<int> hiddenOutputs = {0, 2};
  // program.hideOutputs(hiddenOutputs);

  program.run();

  return 0;
}
