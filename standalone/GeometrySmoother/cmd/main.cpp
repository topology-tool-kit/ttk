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
      = dynamic_cast<ScalarFieldSmoother *>(Program<ttkModule>::ttkModule_);

    if(smoother == nullptr) {
      return 0;
    }

    smoother->setDimensionNumber(3);
    smoother->setInputDataPointer(pointSet_.data());
    smoother->setOutputDataPointer(pointSet_.data());

    // template call based on the triangulation type
    ttkTemplateMacro(triangleMesh_.getType(),
                     (smoother->smooth<float, TTK_TT>(
                       (TTK_TT *)triangleMesh_.getData(), iterationNumber_)));

    return 0;
  }

  int load(const vector<string> &inputPaths) {

    if(inputPaths.empty())
      return -1;
    if(!inputPaths[0].length())
      return -2;

    // choose where to display this message (cout, cerr, a file)
    // choose the priority of this message (1, nearly always displayed,
    // higher values mean lower priorities)
    Debug::printMsg("Reading input mesh...");

    int vertexNumber = 0, triangleNumber = 0;
    string keyword;

    ifstream f(inputPaths[0].data(), ios::in);

    if(!f) {
      Debug::printErr("Cannot open file `" + inputPaths[0] + "'!");
      return -1;
    }

    f >> keyword;

    if(keyword != "OFF") {
      Debug::printErr("Input OFF file `" + inputPaths[0]
                      + "' seems invalid :(");
      return -2;
    }

    f >> vertexNumber;
    f >> triangleNumber;
    f >> keyword;

    pointSet_.resize(3 * vertexNumber);
    triangleSetCo_.resize(3 * triangleNumber);
    triangleSetOff_.resize(triangleNumber + 1);

    for(int i = 0; i < 3 * vertexNumber; i++) {
      f >> pointSet_[i];
    }

    int offId = 0;
    int coId = 0;
    for(int i = 0; i < triangleNumber; i++) {
      int cellSize;
      f >> cellSize;
      if(cellSize != 3) {
        std::cerr << "cell size " << cellSize << " != 3" << std::endl;
        return -3;
      }
      triangleSetOff_[offId++] = coId;
      for(int j = 0; j < 3; j++) {
        int cellId;
        f >> cellId;
        triangleSetCo_[coId++] = cellId;
      }
    }
    triangleSetOff_[offId] = coId; // the last one

    f.close();

    ScalarFieldSmoother *smoother
      = dynamic_cast<ScalarFieldSmoother *>(Program<ttkModule>::ttkModule_);
    triangleMesh_.setInputPoints(vertexNumber, pointSet_.data());
#ifdef TTK_CELL_ARRAY_NEW
    triangleMesh_.setInputCells(
      triangleNumber, triangleSetCo_.data(), triangleSetOff_.data());
#else
    LongSimplexId *triangleSet;
    CellArray::TranslateToFlatLayout(
      triangleSetCo_, triangleSetOff_, triangleSet);
    triangleMesh_.setInputCells(triangleNumber, triangleSet);
#endif

    smoother->preconditionTriangulation(&triangleMesh_);

    Debug::printMsg("done! (read " + std::to_string(vertexNumber)
                    + " vertices, " + std::to_string(triangleNumber)
                    + " triangles)");

    return 0;
  }

  int save() const {

    string fileName(Program<ttkModule>::outputPath_ + ".off");

    ofstream f(fileName.data(), ios::out);

    if(!f) {
      Debug::printErr("Could not write output file `" + fileName + "'!");
      return -1;
    }

    f << "OFF" << endl;
    f << pointSet_.size() / 3 << " " << triangleSetOff_.size() - 1 << " 0"
      << endl;

    for(int i = 0; i < (int)pointSet_.size() / 3; i++) {
      for(int j = 0; j < 3; j++) {
        f << pointSet_[3 * i + j];
        f << " ";
      }
      f << endl;
    }

    for(int i = 0; i < (int)triangleSetOff_.size() - 1; i++) {
      f << "3 ";
      for(int j = 0; j < 3; j++) {
        f << triangleSetCo_[triangleSetOff_[i] + j] << " ";
      }
      f << endl;
    }

    f.close();

    return 0;
  }

protected:
  vector<float> pointSet_;
  vector<long long int> triangleSetCo_;
  vector<long long int> triangleSetOff_;
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
