/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief dummy VTK-free program example that smoothes the geometry of an
/// input surface.

// include the local headers
#include <ProgramBase.h>
#include <ScalarFieldSmoother.h>

using namespace std;
using namespace ttk;

int iterationNumber_ = 1;

template <class ttkModule>
class ExampleProgram : public Program<ttkModule> {

public:
  int execute() {

    ScalarFieldSmoother *smoother
      = (ScalarFieldSmoother *)Program<ttkModule>::ttkModule_;

    smoother->setDimensionNumber(3);
    smoother->setInputDataPointer(pointSet_.data());
    smoother->setOutputDataPointer(pointSet_.data());
    smoother->smooth<float>(iterationNumber_);

    return 0;
  }

  int load(const vector<string> &inputPaths) {

    if(inputPaths.empty())
      return -1;
    if(!inputPaths[0].length())
      return -2;

    {
      stringstream msg;
      msg << "[ExampleProgram] Reading input mesh..." << endl;
      // choose where to display this message (cout, cerr, a file)
      // choose the priority of this message (1, nearly always displayed,
      // higher values mean lower priorities)
      Debug::dMsg(cout, msg.str(), Debug::timeMsg);
    }

    int vertexNumber = 0, triangleNumber = 0;
    string keyword;

    ifstream f(inputPaths[0].data(), ios::in);

    if(!f) {
      stringstream msg;
      msg << "[Editor] Cannot open file `" << inputPaths[0] << "'!" << endl;
      Debug::dMsg(cerr, msg.str(), Debug::fatalMsg);
      return -1;
    }

    f >> keyword;

    if(keyword != "OFF") {
      stringstream msg;
      msg << "[Editor] Input OFF file `" << inputPaths[0]
          << "' seems invalid :(" << endl;
      Debug::dMsg(cerr, msg.str(), Debug::fatalMsg);
      return -2;
    }

    f >> vertexNumber;
    f >> triangleNumber;
    f >> keyword;

    pointSet_.resize(3 * vertexNumber);
    triangleSet_.resize(4 * triangleNumber);

    for(int i = 0; i < 3 * vertexNumber; i++) {
      f >> pointSet_[i];
    }

    for(int i = 0; i < 4 * triangleNumber; i++) {
      f >> triangleSet_[i];
    }

    f.close();

    ScalarFieldSmoother *smoother
      = (ScalarFieldSmoother *)Program<ttkModule>::ttkModule_;

    triangleMesh_.setInputPoints(vertexNumber, pointSet_.data());
    triangleMesh_.setInputCells(triangleNumber, triangleSet_.data());
    smoother->setupTriangulation(&triangleMesh_);

    {
      stringstream msg;
      msg << "[Editor]   done! (read " << vertexNumber << " vertices, "
          << triangleNumber << " triangles)" << endl;
      Debug::dMsg(cout, msg.str(), Debug::timeMsg);
    }

    return 0;
  }

  int save() const {

    string fileName(Program<ttkModule>::outputPath_ + ".off");

    ofstream f(fileName.data(), ios::out);

    if(!f) {
      stringstream msg;
      msg << "[Editor] Could not write output file `" << fileName << "'!"
          << endl;
      Debug::dMsg(cerr, msg.str(), Debug::fatalMsg);
      return -1;
    }

    f << "OFF" << endl;
    f << pointSet_.size() / 3 << " " << triangleSet_.size() / 4 << " 0" << endl;

    for(int i = 0; i < (int)pointSet_.size() / 3; i++) {
      for(int j = 0; j < 3; j++) {
        f << pointSet_[3 * i + j];
        f << " ";
      }
      f << endl;
    }

    for(int i = 0; i < (int)triangleSet_.size() / 4; i++) {
      for(int j = 0; j < 4; j++) {
        f << triangleSet_[4 * i + j];
        f << " ";
      }
      f << endl;
    }

    f.close();

    return 0;
  }

protected:
  vector<float> pointSet_;
  vector<long long int> triangleSet_;
  Triangulation triangleMesh_;
};

int main(int argc, char **argv) {

  ExampleProgram<ScalarFieldSmoother> program;

  // register the arguments to the command line parser
  program.parser_.setArgument(
    "I", &iterationNumber_, "Number of smoothing iterations", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
