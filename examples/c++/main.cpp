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
#include <TopologyToolKit.h>

#include <cassert>
#include <iostream>

int load(const std::string &inputPath,
         std::vector<float> &pointSet,
         std::vector<long long int> &triangleSetCo,
         std::vector<long long int> &triangleSetOff) {

  // load some terrain from some OFF file.

  if(inputPath.empty())
    return -1;

  ttk::Debug dbg;
  dbg.setDebugLevel(ttk::globalDebugLevel_);
  dbg.setDebugMsgPrefix("main::load");
  dbg.printMsg("Reading input mesh...");

  int vertexNumber = 0, triangleNumber = 0;
  std::string keyword;

  std::ifstream f(inputPath.data(), std::ios::in);

  if(!f) {
    dbg.printErr("Cannot read file `" + inputPath + "'!");
    return -1;
  }

  f >> keyword;

  if(keyword != "OFF") {
    dbg.printErr("Input OFF file `" + inputPath + "' seems invalid :(");
    return -2;
  }

  f >> vertexNumber;
  f >> triangleNumber;
  f >> keyword;

  pointSet.resize(3 * vertexNumber);
  triangleSetCo.resize(3 * triangleNumber);
  triangleSetOff.resize(triangleNumber + 1);

  for(int i = 0; i < 3 * vertexNumber; i++) {
    f >> pointSet[i];
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
    triangleSetOff[offId++] = coId;
    for(int j = 0; j < 3; j++) {
      int cellId;
      f >> cellId;
      triangleSetCo[coId++] = cellId;
    }
  }
  triangleSetOff[offId] = coId; // the last one

  f.close();

  dbg.printMsg("... done! (read " + std::to_string(vertexNumber) + " vertices, "
               + std::to_string(triangleNumber) + " triangles)");

  return 0;
}

int save(const std::vector<float> &pointSet,
         const std::vector<long long int> &triangleSetCo,
         const std::vector<long long int> &triangleSetOff,
         const std::string &outputPath) {

  // save the simplified terrain in some OFF file
  std::string fileName(outputPath);

  std::ofstream f(fileName.data(), std::ios::out);

  if(!f) {
    ttk::Debug dbg;
    dbg.setDebugLevel(ttk::globalDebugLevel_);
    dbg.setDebugMsgPrefix("main::save");
    dbg.printErr("Could not write output file `" + fileName + "'!");
    return -1;
  }

  const int nbTriangles = triangleSetOff.size() - 1;

  f << "OFF" << std::endl;
  f << pointSet.size() / 3 << " " << nbTriangles << " 0" << std::endl;

  for(int i = 0; i < (int)pointSet.size() / 3; i++) {
    for(int j = 0; j < 3; j++) {
      f << pointSet[3 * i + j];
      f << " ";
    }
    f << std::endl;
  }

  for(int i = 0; i < nbTriangles; i++) {
    int cellSize = triangleSetOff[i + 1] - triangleSetOff[i];
    assert(cellSize == 3);
    f << cellSize << " ";
    for(int j = triangleSetOff[i]; j < triangleSetOff[i + 1]; j++) {
      f << triangleSetCo[j];
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
  std::vector<long long int> triangleSetCo, triangleSetOff;
  ttk::Triangulation triangulation;

  // load the input
  load(inputFilePath, pointSet, triangleSetCo, triangleSetOff);
  triangulation.setInputPoints(pointSet.size() / 3, pointSet.data());
  long long int triangleNumber = triangleSetOff.size() - 1;
#ifdef TTK_CELL_ARRAY_NEW
  triangulation.setInputCells(
    triangleNumber, triangleSetCo.data(), triangleSetOff.data());
#else
  LongSimplexId *triangleSet;
  CellArray::TranslateToFlatLayout(triangleSetCo, triangleSetOff, triangleSet);
  triangulation.setInputCells(triangleNumber, triangleSet);
#endif

  // NOW, do the TTK processing

  // computing some elevation
  std::vector<float> height(pointSet.size() / 3);
  int vertexId = 0;
  // use the z-coordinate here
  for(int i = 2; i < (int)pointSet.size(); i += 3) {
    height[vertexId] = pointSet[i];
    vertexId++;
  }

  // order array: every vertex sorted according to the elevation field
  std::vector<ttk::SimplexId> order(height.size());
  // precondition/fill in the order array according to the elevation field
  ttk::preconditionOrderArray(height.size(), height.data(), order.data());

  // 2. computing the persistence curve
  ttk::PersistenceCurve curve;
  std::vector<std::pair<float, ttk::SimplexId>> outputCurve;
  curve.preconditionTriangulation(&triangulation);
  curve.setOutputCTPlot(&outputCurve);
  curve.execute<float>(height.data(), order.data(), &triangulation);

  // 3. computing the persitence diagram
  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput;
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(diagramOutput, height.data(), order.data(), &triangulation);

  // 4. selecting the critical point pairs
  std::vector<float> simplifiedHeight = height;
  std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
  for(int i = 0; i < (int)diagramOutput.size(); i++) {
    if(diagramOutput[i].persistence > 0.05) {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(diagramOutput[i].birth);
      authorizedCriticalPoints.push_back(diagramOutput[i].death);
    }
  }

  // 6. simplifying the input data to remove non-persistent pairs
  ttk::TopologicalSimplification simplification;
  simplification.preconditionTriangulation(&triangulation);
  simplification.execute<float>(height.data(), simplifiedHeight.data(),
                                authorizedCriticalPoints.data(), order.data(),
                                simplifiedOrder.data(),
                                authorizedCriticalPoints.size(), triangulation);

  // assign the simplified values to the input mesh
  for(int i = 0; i < (int)simplifiedHeight.size(); i++) {
    pointSet[3 * i + 2] = simplifiedHeight[i];
  }

  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  ttk::MorseSmaleComplex::OutputCriticalPoints outCriticalPoints{};
  // 1-separatrices
  ttk::MorseSmaleComplex::Output1Separatrices out1Separatrices{};
  // 2-separatrices
  ttk::MorseSmaleComplex::Output2Separatrices out2Separatrices{};
  // segmentation
  std::vector<ttk::SimplexId> ascendingSegmentation(
    triangulation.getNumberOfVertices(), -1),
    descendingSegmentation(triangulation.getNumberOfVertices(), -1),
    mscSegmentation(triangulation.getNumberOfVertices(), -1);
  ttk::MorseSmaleComplex::OutputManifold outSegmentation{
    ascendingSegmentation.data(), descendingSegmentation.data(),
    mscSegmentation.data()};

  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.execute(
    outCriticalPoints, out1Separatrices, out2Separatrices, outSegmentation,
    simplifiedHeight.data(), simplifiedOrder.data(), triangulation);

  // save the output
  save(pointSet, triangleSetCo, triangleSetOff, "output.off");

  return 0;
}
