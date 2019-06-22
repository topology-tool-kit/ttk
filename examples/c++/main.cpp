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
#include <MorseSmaleComplex.h>
#include <PersistenceCurve.h>
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>

int load(const std::string &inputPath,
         std::vector<float> &pointSet,
         std::vector<long long int> &triangleSet) {

  // load some terrain from some OFF file.

  if(inputPath.empty())
    return -1;

  ttk::Debug d;

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
  std::vector<int> offsets(height.size());
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
  curve.execute<float, int>();

  // 3. computing the persitence diagram
  ttk::PersistenceDiagram diagram;
  std::vector<std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,
                         ttk::CriticalType, float, ttk::SimplexId>>
    diagramOutput;
  diagram.setupTriangulation(&triangulation);
  diagram.setInputScalars(height.data());
  diagram.setInputOffsets(offsets.data());
  diagram.setOutputCTDiagram(&diagramOutput);
  diagram.execute<float, int>();

  // 4. selecting the critical point pairs
  std::vector<float> simplifiedHeight = height;
  std::vector<int> authorizedCriticalPoints, simplifiedOffsets = offsets;
  for(int i = 0; i < (int)diagramOutput.size(); i++) {
    double persistence = std::get<4>(diagramOutput[i]);
    if(persistence > 0.05) {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(std::get<0>(diagramOutput[i]));
      authorizedCriticalPoints.push_back(std::get<2>(diagramOutput[i]));
    }
  }

  // 6. simplifying the input data to remove non-persistent pairs
  ttk::TopologicalSimplification simplification;
  simplification.setupTriangulation(&triangulation);
  simplification.setInputScalarFieldPointer(height.data());
  simplification.setInputOffsetScalarFieldPointer(offsets.data());
  simplification.setOutputOffsetScalarFieldPointer(simplifiedOffsets.data());
  simplification.setOutputScalarFieldPointer(simplifiedHeight.data());
  simplification.setConstraintNumber(authorizedCriticalPoints.size());
  simplification.setVertexIdentifierScalarFieldPointer(
    authorizedCriticalPoints.data());
  simplification.execute<float, int>();

  // assign the simplified values to the input mesh
  for(int i = 0; i < (int)simplifiedHeight.size(); i++) {
    pointSet[3 * i + 2] = simplifiedHeight[i];
  }

  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  int criticalPoints_numberOfPoints{};
  std::vector<float> criticalPoints_points;
  std::vector<char> criticalPoints_points_cellDimensions;
  std::vector<int> criticalPoints_points_cellIds;
  std::vector<char> criticalPoints_points_isOnBoundary;
  std::vector<float> criticalPoints_points_cellScalars;
  std::vector<int> criticalPoints_points_PLVertexIdentifiers;
  std::vector<int> criticalPoints_points_manifoldSize;
  // 1-separatrices
  int separatrices1_numberOfPoints{};
  std::vector<float> separatrices1_points;
  std::vector<char> separatrices1_points_smoothingMask;
  std::vector<char> separatrices1_points_cellDimensions;
  std::vector<int> separatrices1_points_cellIds;
  int separatrices1_numberOfCells{};
  std::vector<int> separatrices1_cells;
  std::vector<int> separatrices1_cells_sourceIds;
  std::vector<int> separatrices1_cells_destinationIds;
  std::vector<int> separatrices1_cells_separatrixIds;
  std::vector<char> separatrices1_cells_separatrixTypes;
  std::vector<char> separatrices1_cells_isOnBoundary;
  std::vector<float> separatrices1_cells_separatrixFunctionMaxima;
  std::vector<float> separatrices1_cells_separatrixFunctionMinima;
  std::vector<float> separatrices1_cells_separatrixFunctionDiffs;
  // segmentation
  std::vector<int> ascendingSegmentation(
    triangulation.getNumberOfVertices(), -1),
    descendingSegmentation(triangulation.getNumberOfVertices(), -1),
    mscSegmentation(triangulation.getNumberOfVertices(), -1);
  morseSmaleComplex.setupTriangulation(&triangulation);
  morseSmaleComplex.setInputScalarField(simplifiedHeight.data());
  morseSmaleComplex.setInputOffsets(simplifiedOffsets.data());
  morseSmaleComplex.setOutputMorseComplexes(ascendingSegmentation.data(),
                                            descendingSegmentation.data(),
                                            mscSegmentation.data());
  morseSmaleComplex.setOutputCriticalPoints(
    &criticalPoints_numberOfPoints, &criticalPoints_points,
    &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
    &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
    &criticalPoints_points_PLVertexIdentifiers,
    &criticalPoints_points_manifoldSize);
  morseSmaleComplex.setOutputSeparatrices1(
    &separatrices1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &separatrices1_numberOfCells,
    &separatrices1_cells, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes,
    &separatrices1_cells_separatrixFunctionMaxima,
    &separatrices1_cells_separatrixFunctionMinima,
    &separatrices1_cells_separatrixFunctionDiffs,
    &separatrices1_cells_isOnBoundary);

  morseSmaleComplex.execute<float, int>();

  // save the output
  save(pointSet, triangleSet, "output.off");

  return 0;
}
