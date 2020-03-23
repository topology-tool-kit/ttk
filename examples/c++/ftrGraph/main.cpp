/// \ingroup examples
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2017.
///
/// \brief Minimalist C++-only TTK example pipeline, including:
///  -# The computation of a persistence curve
///  -# The computation of a persistence diagram
///  -# The selection of the most persistent pairs of the diagram
///  -# The pre-simplification of the data according to this selection
///  -# The computation of the Morse-Smale complex on this simplified data
///  -# The storage of the output of this pipeline to disk.
///
/// This example reproduces the Figure 1 of the TTK companion paper:
/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
/// of IEEE VIS 2017.
///
/// See the individual VTK wrappers (core/vtk/) to see how to use each ttk::base
/// (C++-only) TTK component.

// include the local headers
#include <CommandLineParser.h>
#include <PersistenceCurve.h>
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>
#include <FTRGraph.h>

int load(const std::string &inputPath,
         std::vector<float> &pointSet,
         std::vector<long long int> &triangleSet) {

  // load some terrain from some OFF file.

  if(inputPath.empty())
    return -1;

  ttk::Debug d;
  d.setDebugLevel(ttk::globalDebugLevel_);

  {
    std::stringstream msg;
    msg << "[main::load] Reading input mesh..." << std::endl;
    // choose where to display this message (std::cout, std::cerr, a file)
    // choose the priority of this message (1, nearly always displayed,
    // higher values mean lower priorities)
    d.dMsg(std::cout, msg.str(), d.timeMsg);
  }

  int vertexNumber = 0, triangleNumber = 0;
  std::string keyword;

  std::ifstream f(inputPath.data(), std::ios::in);

  if(!f) {
    std::stringstream msg;
    msg << "[main::load] Cannot open file `" << inputPath << "'!" << std::endl;
    d.dMsg(std::cerr, msg.str(), d.fatalMsg);
    return -1;
  }

  f >> keyword;

  if(keyword != "OFF") {
    std::stringstream msg;
    msg << "[main::load] Input OFF file `" << inputPath << "' seems invalid :("
        << std::endl;
    d.dMsg(std::cerr, msg.str(), d.fatalMsg);
    return -2;
  }

  f >> vertexNumber;
  f >> triangleNumber;
  f >> keyword;

  pointSet.resize(3 * vertexNumber);
  triangleSet.resize(4 * triangleNumber);

  for(int i = 0; i < 3 * vertexNumber; i++) {
    f >> pointSet[i];
  }

  for(int i = 0; i < 4 * triangleNumber; i++) {
    f >> triangleSet[i];
  }

  f.close();

  {
    std::stringstream msg;
    msg << "[main::load]   done! (read " << vertexNumber << " vertices, "
        << triangleNumber << " triangles)" << std::endl;
    d.dMsg(std::cout, msg.str(), d.timeMsg);
  }

  return 0;
}

int save(const std::vector<float> &pointSet,
         const std::vector<long long int> &triangleSet,
         const std::string &outputPath) {

  // save the simplified terrain in some OFF file
  ttk::Debug d;
  d.setDebugLevel(ttk::globalDebugLevel_);

  std::string fileName(outputPath);

  std::ofstream f(fileName.data(), std::ios::out);

  if(!f) {
    std::stringstream msg;
    msg << "[main::save] Could not write output file `" << fileName << "'!"
        << std::endl;
    d.dMsg(std::cerr, msg.str(), d.fatalMsg);
    return -1;
  }

  f << "OFF" << std::endl;
  f << pointSet.size() / 3 << " " << triangleSet.size() / 4 << " 0"
    << std::endl;

  for(int i = 0; i < (int)pointSet.size() / 3; i++) {
    for(int j = 0; j < 3; j++) {
      f << pointSet[3 * i + j];
      f << " ";
    }
    f << std::endl;
  }

  for(int i = 0; i < (int)triangleSet.size() / 4; i++) {
    for(int j = 0; j < 4; j++) {
      f << triangleSet[4 * i + j];
      f << " ";
    }
    f << std::endl;
  }

  f.close();

  return 0;
}

int main(int argc, char **argv) {

  std::string inputFilePath;
  ttk::CommandLineParser parser;

  ttk::globalDebugLevel_ = 3;

  // register the arguments to the command line parser
  parser.setArgument("i", &inputFilePath, "Path to input OFF file");
  // parse
  parser.parse(argc, argv);

  std::vector<float> pointSet;
  std::vector<long long int> triangleSet;
  ttk::Triangulation triangulation;

  // load the input
  load(inputFilePath, pointSet, triangleSet);
  triangulation.setInputPoints(pointSet.size() / 3, pointSet.data());
  triangulation.setInputCells(triangleSet.size() / 4, triangleSet.data());

  // NOW, do the TTK processing

  // computing some elevation
  std::vector<float> height(pointSet.size() / 3);
  std::vector<ttk::SimplexId> offsets(height.size());
  int vertexId = 0;
  // use the z-coordinate here
  for(int i = 2; i < (int)pointSet.size(); i += 3) {
    height[vertexId] = pointSet[i];
    offsets[vertexId] = vertexId;
    vertexId++;
  }

  // 2. computing the persistence curve
  ttk::PersistenceCurve curve;
  std::vector<std::pair<float, ttk::SimplexId>> outputCurve;
  curve.setupTriangulation(&triangulation);
  curve.setInputScalars(height.data());
  curve.setInputOffsets(offsets.data());
  curve.setOutputCTPlot(&outputCurve);
  curve.execute<float, ttk::SimplexId>();

  // 3. computing the persitence diagram
  ttk::PersistenceDiagram diagram;
  std::vector<std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,
                         ttk::CriticalType, float, ttk::SimplexId>>
    diagramOutput;
  diagram.setupTriangulation(&triangulation);
  diagram.setInputScalars(height.data());
  diagram.setInputOffsets(offsets.data());
  diagram.setOutputCTDiagram(&diagramOutput);
  diagram.execute<float, ttk::SimplexId>();

  // 4. selecting the critical point pairs
  std::vector<float> simplifiedHeight = height;
  std::vector<ttk::SimplexId> authorizedCriticalPoints,
    simplifiedOffsets = offsets;
  for(int i = 0; i < (int)diagramOutput.size(); i++) {
    double persistence = std::get<4>(diagramOutput[i]);
    if(persistence > 0.05) {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(std::get<0>(diagramOutput[i]));
      authorizedCriticalPoints.push_back(std::get<2>(diagramOutput[i]));
    }
  }

    // 6. compute the Reeb Graph
    ttk::ftr::FTRGraph<float> ftrGraph;
    ftrGraph.setDebugLevel(3);
    ftrGraph.setThreadNumber(1);
    ftrGraph.setupTriangulation(&triangulation);
    ftrGraph.setScalars(simplifiedHeight.data());
    ftrGraph.setVertexSoSoffsets(&simplifiedOffsets);
    ttk::ftr::Params graphParams;
    graphParams.singleSweep = false;
    graphParams.segm = true;
    graphParams.normalize = true;
    graphParams.advStats = true;
    graphParams.samplingLvl = 0;
    graphParams.threadNumber = 1;
    graphParams.debugLevel = 3;
    ftrGraph.setParams(graphParams);
    ftrGraph.build();

  // save the output
  save(pointSet, triangleSet, "output.off");

  return 0;
}
