#include <MandatoryCriticalPoints.h>

using namespace ttk;

MandatoryCriticalPoints::MandatoryCriticalPoints() {
  this->setDebugMsgPrefix("MandatoryCriticalPoints");
  upperJoinTree_.setDebugLevel(debugLevel_);
  lowerJoinTree_.setDebugLevel(debugLevel_);
  upperSplitTree_.setDebugLevel(debugLevel_);
  lowerSplitTree_.setDebugLevel(debugLevel_);
}

void MandatoryCriticalPoints::flush() {
  inputUpperBoundField_ = NULL;
  inputLowerBoundField_ = NULL;
  outputMandatoryMinimum_ = NULL;
  outputMandatoryJoinSaddle_ = NULL;
  outputMandatorySplitSaddle_ = NULL;
  outputMandatoryMaximum_ = NULL;
  vertexNumber_ = 0;
  upperJoinTree_.flush();
  lowerJoinTree_.flush();
  upperSplitTree_.flush();
  lowerSplitTree_.flush();
  normalizedThreshold_ = 0.0;
  vertexPositions_.clear();
  vertexSoSoffsets_.clear();
  upperVertexScalars_.clear();
  lowerVertexScalars_.clear();
  upperMinimumList_.clear();
  lowerMinimumList_.clear();
  upperMaximumList_.clear();
  lowerMaximumList_.clear();
  mandatoryMinimumVertex_.clear();
  mandatoryMaximumVertex_.clear();
  mandatoryMinimumInterval_.clear();
  mandatoryMaximumInterval_.clear();
  mandatoryJoinSaddleVertex_.clear();
  mandatorySplitSaddleVertex_.clear();
  mergedMaximaId_.clear();
  mergedMinimaId_.clear();
  mdtMinJoinSaddlePair_.clear();
  mdtMaxSplitSaddlePair_.clear();
  isMdtMinimumSimplified_.clear();
  isMdtJoinSaddleSimplified_.clear();
  isMdtSplitSaddleSimplified_.clear();
  isMdtMaximumSimplified_.clear();
  mdtMinimumParentSaddleId_.clear();
  mdtJoinSaddleParentSaddleId_.clear();
  mdtSplitSaddleParentSaddleId_.clear();
  mdtMaximumParentSaddleId_.clear();
  mdtJoinTree_.clear();
  mdtSplitTree_.clear();
  mdtJoinTreePointComponentId_.clear();
  mdtSplitTreePointComponentId_.clear();
  mdtJoinTreePointType_.clear();
  mdtSplitTreePointType_.clear();
  mdtJoinTreePointLowInterval_.clear();
  mdtSplitTreePointLowInterval_.clear();
  mdtJoinTreePointUpInterval_.clear();
  mdtSplitTreePointUpInterval_.clear();
  mdtJoinTreePointXCoord_.clear();
  mdtSplitTreePointXCoord_.clear();
  mdtJoinTreePointYCoord_.clear();
  mdtSplitTreePointYCoord_.clear();
  mandatoryMaximumComponentVertices_.clear();
  mandatoryMinimumComponentVertices_.clear();
  mandatoryJoinSaddleComponentVertices_.clear();
  mandatorySplitSaddleComponentVertices_.clear();
}

int MandatoryCriticalPoints::buildMandatoryTree(
  const TreeType treeType,
  Graph &mdtTree,
  std::vector<int> &mdtTreePointComponentId,
  std::vector<PointType> &mdtTreePointType,
  std::vector<double> &mdtTreePointLowInterval,
  std::vector<double> &mdtTreePointUpInterval,
  std::vector<int> &mdtTreeEdgeSwitchable,
  const std::vector<int> &mdtExtremumParentSaddle,
  const std::vector<int> &mdtSaddleParentSaddle,
  const std::vector<bool> &isExtremumSimplified,
  const std::vector<bool> &isSaddleSimplified,
  const std::vector<std::pair<double, double>> &extremumInterval,
  const std::vector<std::pair<int, int>> &mandatorySaddleVertices,
  const int extremaNumber,
  const int saddleNumber,
  const PointType extremumType,
  const PointType saddleType,
  const PointType otherExtremumType,
  const double globalOtherExtremumValue) const {

  const std::string treeName
    = treeType == TreeType::JoinTree ? "join tree" : "split tree";

  /* Preliminaries operations */
  mdtTree.clear();
  mdtTreePointComponentId.clear();
  mdtTreePointComponentId.reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointType.clear();
  mdtTreePointType.reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointLowInterval.clear();
  mdtTreePointLowInterval.reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointUpInterval.clear();
  mdtTreePointUpInterval.reserve(extremaNumber + saddleNumber + 1);
  mdtTreeEdgeSwitchable.clear();

  /* Graph Object Building */
  std::vector<int> saddleGraphVertex(saddleNumber, -1);
  // Create and connect all extrema to their saddle
  for(int i = 0; i < extremaNumber; i++) {
    // If simplified, do nothing and go to the next extremum
    if(isExtremumSimplified[i])
      continue;
    // New point in the graph
    int extremumGraphPoint = mdtTree.addVertex();
    mdtTreePointComponentId.push_back(i);
    mdtTreePointType.push_back(extremumType);

    mdtTreePointLowInterval.push_back(extremumInterval[i].first);
    mdtTreePointUpInterval.push_back(extremumInterval[i].second);
    // Look for the saddle and connect if there is one
    int parentSaddle = mdtExtremumParentSaddle[i];
    // If no parent saddle, end the loop
    if(parentSaddle == -1)
      break;
    // Create the saddle (if not already)
    if(saddleGraphVertex[parentSaddle] == -1) {
      saddleGraphVertex[parentSaddle] = mdtTree.addVertex();
      mdtTreePointComponentId.push_back(parentSaddle);
      mdtTreePointType.push_back(saddleType);
      int lowerVertex = mandatorySaddleVertices[parentSaddle].first;
      int upperVertex = mandatorySaddleVertices[parentSaddle].second;
      mdtTreePointLowInterval.push_back(lowerVertexScalars_[lowerVertex]);
      mdtTreePointUpInterval.push_back(upperVertexScalars_[upperVertex]);
    }
    // Connect the extrema and saddle
    mdtTree.addEdge(saddleGraphVertex[parentSaddle], extremumGraphPoint);
    // Not switchable (saddle -> extremum)
    mdtTreeEdgeSwitchable.push_back(0);
  }

  // If there is no extremum (=> no saddles), do nothing (empty graph)
  if(mdtTree.getNumberOfVertices() == 0) {
    this->printWrn("Mandatory " + treeName + " empty.");
    return 0;
  }

  // If there is only one extremum (=> no saddles), connect the extremum to the
  // other global extremum
  if(mdtTree.getNumberOfVertices() == 1) {
    mdtTree.addVertex();
    mdtTreePointComponentId.push_back(-1);
    mdtTreePointType.push_back(otherExtremumType);
    mdtTreePointLowInterval.push_back(globalOtherExtremumValue);
    mdtTreePointUpInterval.push_back(globalOtherExtremumValue);
    mdtTree.addEdge(1, 0);
    mdtTreeEdgeSwitchable.push_back(0);
    this->printWrn("Mandatory " + treeName + " with only one minimum.");

  } else {
    // Continue to add the remaining saddles
    // Create and connect all remaining saddles
    for(int i = 0; i < saddleNumber; i++) {

      // If simplified, do nothing and go to the next saddle
      if(isSaddleSimplified[i]) {
        continue;
      }
      // Create the graph point if not already
      if(saddleGraphVertex[i] == -1) {
        saddleGraphVertex[i] = mdtTree.addVertex();
        mdtTreePointComponentId.push_back(i);
        mdtTreePointType.push_back(saddleType);
        int lowerVertex = mandatorySaddleVertices[i].first;
        int upperVertex = mandatorySaddleVertices[i].second;
        mdtTreePointLowInterval.push_back(lowerVertexScalars_[lowerVertex]);
        mdtTreePointUpInterval.push_back(upperVertexScalars_[upperVertex]);
      }
      // Look for the saddle above and connect if there is one
      int parentSaddle = mdtSaddleParentSaddle[i];
      // If the parent is different from the saddle itself, create it and
      // connect
      if(parentSaddle != i) {
        // Create the graph point if not already
        if(saddleGraphVertex[parentSaddle] == -1) {
          saddleGraphVertex[parentSaddle] = mdtTree.addVertex();
          mdtTreePointComponentId.push_back(parentSaddle);
          mdtTreePointType.push_back(saddleType);
          int lowerVertex = mandatorySaddleVertices[parentSaddle].first;
          int upperVertex = mandatorySaddleVertices[parentSaddle].second;
          mdtTreePointLowInterval.push_back(lowerVertexScalars_[lowerVertex]);
          mdtTreePointUpInterval.push_back(upperVertexScalars_[upperVertex]);
        }
        // Connect the two saddles
        mdtTree.addEdge(saddleGraphVertex[parentSaddle], saddleGraphVertex[i]);
        // Test if switchable : parentSaddle and i
        mdtTreeEdgeSwitchable.push_back(
          areSaddlesSwitchables(treeType, parentSaddle, i));
      } else { // It is the root, create the global extremum to connect with it
        int globalOtherExtremumGraphPoint = mdtTree.addVertex();
        mdtTreePointComponentId.push_back(-1);
        mdtTreePointType.push_back(otherExtremumType);
        mdtTreePointLowInterval.push_back(globalOtherExtremumValue);
        mdtTreePointUpInterval.push_back(globalOtherExtremumValue);
        mdtTree.addEdge(globalOtherExtremumGraphPoint, saddleGraphVertex[i]);
        mdtTreeEdgeSwitchable.push_back(0);
      }
    }
  }

  /* Debug Messages */
  this->printMsg("Building of the mandatory " + treeName + "");
  this->printMsg({{"#Points", std::to_string(mdtTree.getNumberOfVertices())},
                  {"#Edges", std::to_string(mdtTree.getNumberOfEdges())}});

  this->printMsg("List of mandatory " + treeName + " graph points:",
                 debug::Priority::DETAIL);
  for(int i = 0; i < mdtTree.getNumberOfVertices(); i++) {
    std::stringstream msg{};
    msg << "  (" << std::setw(3) << std::right << i << ") ";
    msg << std::setw(12) << std::right;
    switch(mdtTreePointType[i]) {
      case PointType::Minimum:
        msg << "Minimum";
        break;
      case PointType::Maximum:
        msg << "Maximum";
        break;
      case PointType::JoinSaddle:
        msg << "Join Saddle";
        break;
      case PointType::SplitSaddle:
        msg << "Split Saddle";
        break;
      default:
        break;
    }
    msg << "  id = " << std::setw(3) << mdtTreePointComponentId[i];
    msg << "  I = [ " << mdtTreePointLowInterval[i];
    msg << " ; " << mdtTreePointUpInterval[i];
    msg << " ]";
    this->printMsg(msg.str(), debug::Priority::DETAIL);
  }

  this->printMsg(
    "List of mandatory " + treeName + " graph edges:", debug::Priority::DETAIL);
  for(int i = 0; i < mdtTree.getNumberOfEdges(); i++) {
    std::stringstream msg{};
    msg << "  (" << std::setw(3) << std::right << i << ") ";
    msg << "Points  ";
    msg << std::setw(3) << std::right
        << mdtTree.getEdge(i).getVertexIdx().first;
    msg << " -> ";
    msg << std::setw(3) << std::left
        << mdtTree.getEdge(i).getVertexIdx().second;
    this->printMsg(msg.str(), debug::Priority::DETAIL);
  }

  return 0;
}

int MandatoryCriticalPoints::buildPairs(
  const TreeType treeType,
  const std::vector<std::pair<int, int>> &saddleList,
  const std::vector<std::vector<int>> &mergedExtrema,
  const std::vector<std::pair<double, double>> &extremumInterval,
  SubLevelSetTree &lowerTree,
  SubLevelSetTree &upperTree,
  std::vector<std::pair<std::pair<int, int>, double>> &extremaSaddlePair)
  const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(lowerTree.isJoinTree() != upperTree.isJoinTree())
    return -1;
#endif

  // List of pairs
  // .first.first = saddle id
  // .first.second = extremum id
  // .second = metric d(M,S)
  extremaSaddlePair.clear();
  // Build list of pairs
  for(size_t i = 0; i < mergedExtrema.size(); i++) {
    for(size_t j = 0; j < mergedExtrema[i].size(); j++) {
      // Build a pair (Si,Mj)
      std::pair<std::pair<int, int>, double> constructedPair{
        {mergedExtrema[i][j], i}, 0.0};
      // Get the values to compute d(Mj,Si)
      double lowerValue = 0;
      double upperValue = 0;
      if(lowerTree.isJoinTree()) {
        lowerValue = extremumInterval[constructedPair.first.first].first;
        // lowerTree->getVertexScalar(extremumList[(*mergedExtrema)[i][j]],
        // lowerValue);
        upperTree.getVertexScalar(saddleList[i].second, upperValue);
      } else {
        lowerTree.getVertexScalar(saddleList[i].first, lowerValue);
        // upperTree->getVertexScalar(extremumList[(*mergedExtrema)[i][j]],
        // upperValue);
        upperValue = extremumInterval[constructedPair.first.first].second;
      }
      // Evaluate d(Si,Mj)
      constructedPair.second = fabs(upperValue - lowerValue);
      // Add the pair to the list
      extremaSaddlePair.push_back(constructedPair);
    }
  }
  // Sort pair by increasing value of d(S,M) (.second)
  criticalPointPairComparaison pairComparaison;
  std::sort(
    extremaSaddlePair.begin(), extremaSaddlePair.end(), pairComparaison);

  const std::string treeName = treeType == TreeType::JoinTree
                                 ? "minimum - join saddle"
                                 : "maximum - split saddle";
  this->printMsg(
    "List of mandatory " + treeName + " pairs:", debug::Priority::DETAIL);
  for(size_t i = 0; i < extremaSaddlePair.size(); i++) {
    std::stringstream msg;
    msg << "  (" << std::setw(3) << extremaSaddlePair[i].first.first << ";";
    msg << std::setw(3) << extremaSaddlePair[i].first.second << ")";
    msg << " -> d = " << extremaSaddlePair[i].second;
    this->printMsg(msg.str(), debug::Priority::DETAIL);
  }

  return 0;
}

int MandatoryCriticalPoints::computePlanarLayout(
  const TreeType &treeType,
  const Graph &mdtTree,
  const std::vector<PointType> &mdtTreePointType,
  const std::vector<double> &mdtTreePointLowInterval,
  const std::vector<double> &mdtTreePointUpInterval,
  std::vector<double> &xCoord,
  std::vector<double> &yCoord) const {

  /* Informations */
  const int numberOfPoints = mdtTree.getNumberOfVertices();
  const double rangeMin = getGlobalMinimum();
  const double rangeMax = getGlobalMaximum();
  const double range = rangeMax - rangeMin;

  // Get the root
  int rootGraphPointId = -1;
  if(treeType == TreeType::JoinTree) {
    for(size_t i = 0; i < mdtTreePointType.size(); i++) {
      if(mdtTreePointType[i] == PointType::Maximum) {
        rootGraphPointId = i;
        break;
      }
    }
  } else {
    for(size_t i = 0; i < mdtTreePointType.size(); i++) {
      if(mdtTreePointType[i] == PointType::Minimum) {
        rootGraphPointId = i;
        break;
      }
    }
  }
  if(rootGraphPointId == -1) {
    return -1;
  }

  // Root down point
  int downRootPointId = -1;
  if(mdtTree.getVertex(rootGraphPointId).getNumberOfEdges() == 1) {
    int edgeId = mdtTree.getVertex(rootGraphPointId).getEdgeIdx(0);
    if(mdtTree.getEdge(edgeId).getVertexIdx().first == rootGraphPointId) {
      downRootPointId = mdtTree.getEdge(edgeId).getVertexIdx().second;
    }
  }
  if(downRootPointId == -1) {
    return -2;
  }

  // Resize of outputs
  xCoord.resize(numberOfPoints);
  yCoord.resize(numberOfPoints);

  // Graph Point x intervals
  std::vector<std::pair<double, double>> xInterval(numberOfPoints);
  xInterval[rootGraphPointId].first = -0.5 * (double)numberOfPoints;
  xInterval[rootGraphPointId].second = 0.5 * (double)numberOfPoints;
  // Order
  std::vector<std::pair<int, double>> xOrder(numberOfPoints);

  std::queue<int> graphPointQueue;
  graphPointQueue.push(rootGraphPointId);

  while(!graphPointQueue.empty()) {

    int graphPoint = graphPointQueue.front();
    graphPointQueue.pop();

    // Number Of Down Edges
    int numberOfEdges = mdtTree.getVertex(graphPoint).getNumberOfEdges();
    int numberOfDownEdges = 0;
    for(int i = 0; i < numberOfEdges; i++) {
      int edgeId = mdtTree.getVertex(graphPoint).getEdgeIdx(i);
      if(mdtTree.getEdge(edgeId).getVertexIdx().first == graphPoint) {
        numberOfDownEdges++;
      }
    }
    // X Order
    xOrder[graphPoint].first = graphPoint;
    xOrder[graphPoint].second
      = (xInterval[graphPoint].second + xInterval[graphPoint].first) * 0.5;
    // Continue to next point in queue if it's a leaf
    if(numberOfDownEdges == 0)
      continue;
    // Down node count
    int downPointCount = 0;
    // Interval splitting for down points
    for(int i = 0; i < numberOfEdges; i++) {
      int edgeId = mdtTree.getVertex(graphPoint).getEdgeIdx(i);
      int downGraphPoint = -1;
      if(mdtTree.getEdge(edgeId).getVertexIdx().first == graphPoint) {
        downGraphPoint = mdtTree.getEdge(edgeId).getVertexIdx().second;
      }
      // Next edge if it is not a down edge
      if(downGraphPoint == -1)
        continue;
      // Interval for the down graph
      xInterval[downGraphPoint].first
        = xInterval[graphPoint].first
          + ((xInterval[graphPoint].second - xInterval[graphPoint].first)
             * (double)downPointCount / (double)numberOfDownEdges);
      xInterval[downGraphPoint].second
        = xInterval[graphPoint].first
          + ((xInterval[graphPoint].second - xInterval[graphPoint].first)
             * (double)(downPointCount + 1) / (double)numberOfDownEdges);
      // Count the point and add it to the queue
      downPointCount++;
      graphPointQueue.push(downGraphPoint);
    }
  }

  /* Y coordinates */
  for(int i = 0; i < numberOfPoints; i++) {
    yCoord[i]
      = ((0.5 * (mdtTreePointLowInterval[i] + mdtTreePointUpInterval[i]))
         - rangeMin)
        / range;
  }

  /* X coordinates */
  // Sorting
  pairComparaison xCoordCmp;
  std::sort(xOrder.begin(), xOrder.end(), xCoordCmp);
  // X coordinates for all points except root (global extremum)
  int pointCount = 0;
  for(int i = 0; i < numberOfPoints; i++) {
    if(xOrder[i].first != rootGraphPointId) {
      xCoord[xOrder[i].first]
        = (double)pointCount * (1.0 / ((int)numberOfPoints - 2.0));
      pointCount++;
    }
  }
  // X coordinate for the root
  xCoord[rootGraphPointId] = xCoord[downRootPointId];

  const std::string treeName
    = treeType == TreeType::JoinTree ? "join tree" : "split tree";

  this->printMsg("Planar layout for mandatory " + treeName + ": root point ("
                   + std::to_string(downRootPointId) + ")",
                 debug::Priority::DETAIL);

  for(int i = 0; i < numberOfPoints; i++) {
    std::stringstream msg;
    msg << "  (" << i << ") :"
        << "  x = " << xCoord[i] << "  y = " << yCoord[i];
    this->printMsg(msg.str(), debug::Priority::DETAIL);
  }

  return 0;
}

int MandatoryCriticalPoints::computeExtremumComponent(
  const PointType &pointType,
  const SubLevelSetTree &tree,
  const int seedVertexId,
  const std::vector<double> &vertexScalars,
  std::vector<int> &componentVertexList) const {

  const double value = vertexScalars[seedVertexId];

  // Clear the list
  componentVertexList.clear();
  // Get the super arc id of the vertex
  int superArcId = getVertexSuperArcId(seedVertexId, &tree);
  // Get the sub tree root super arc
  int rootSuperArcId = getSubTreeRootSuperArcId(&tree, superArcId, value);
  // Get the list of the sub tree super arc ids
  std::vector<int> subTreeSuperArcId;
  getSubTreeSuperArcIds(&tree, rootSuperArcId, subTreeSuperArcId);
  // Comparaison
  int sign = (pointType == PointType::Minimum) ? 1 : -1;
  // Compute each super arc
  for(int i = 0; i < (int)subTreeSuperArcId.size(); i++) {
    const SuperArc *superArc = tree.getSuperArc(subTreeSuperArcId[i]);
    const int numberOfRegularNodes = superArc->getNumberOfRegularNodes();
    // Root super arc and others treated differently
    if(subTreeSuperArcId[i] == rootSuperArcId) {
      // Test the value for each regular node
      for(int j = 0; j < numberOfRegularNodes; j++) {
        int regularNodeId = superArc->getRegularNodeId(j);
        double nodeScalar = tree.getNodeScalar(regularNodeId);
        if(!((sign * nodeScalar) > (sign * value))) {
          int vertexId = tree.getNode(regularNodeId)->getVertexId();
          componentVertexList.push_back(vertexId);
        }
      }
      // Down node
      int downNodeId = superArc->getDownNodeId();
      double nodeScalar = tree.getNodeScalar(downNodeId);
      if(!((sign * nodeScalar) > (sign * value))) {
        int vertexId = tree.getNode(downNodeId)->getVertexId();
        componentVertexList.push_back(vertexId);
      }
    } else {
      // Take all regular nodes
      for(int j = 0; j < numberOfRegularNodes; j++) {
        int regularNodeId = superArc->getRegularNodeId(j);
        int vertexId = tree.getNode(regularNodeId)->getVertexId();
        componentVertexList.push_back(vertexId);
      }
      // Take down node
      int downNodeId = superArc->getDownNodeId();
      int vertexId = tree.getNode(downNodeId)->getVertexId();
      componentVertexList.push_back(vertexId);
    }
  }
  return 0;
}

int MandatoryCriticalPoints::enumerateMandatoryExtrema(
  const PointType pointType,
  SubLevelSetTree &firstTree,
  SubLevelSetTree &secondTree,
  std::vector<int> &mandatoryExtremum,
  std::vector<std::pair<double, double>> &criticalInterval) const {

  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!(firstTree.isJoinTree() != firstTree.isSplitTree()))
    return -1;
  if(!(secondTree.isJoinTree() != secondTree.isSplitTree()))
    return -2;
  if(!((firstTree.isJoinTree() && secondTree.isJoinTree())
       || (firstTree.isSplitTree() && secondTree.isSplitTree())))
    return -3;
#endif

  // Extremum list in the first tree
  const std::vector<int> &extremumList = *firstTree.getExtremumList();
  // Some tmp variables
  size_t extremumNumber = extremumList.size();
  std::vector<int> vertexId(extremumNumber);
  std::vector<double> vertexValue(extremumNumber);
  std::vector<int> superArcId(extremumNumber);
  // vector<int> rootSuperArcId(extremumNumber);
  std::vector<int> subTreeSuperArcId; // vector<vector<int> >
                                      // subTreeSuperArcId(extremumNumber);
  // Clear the list of mandatory extrema
  mandatoryExtremum.clear();
  criticalInterval.clear();
  // To mark the super arc in the second tree as already used for a mandatory
  // critical component
  std::vector<bool> isSuperArcAlreadyVisited(
    secondTree.getNumberOfSuperArcs(), false);

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  for(size_t i = 0; i < extremumNumber; i++) {
    // Mandatory until proven otherwise
    bool isMandatory = true;
    // Vertex Id (first tree)
    vertexId[i] = extremumList[i];
    // Vertex Value (first tree)
    firstTree.getVertexScalar(vertexId[i], vertexValue[i]);
    // Super Arc Id (second tree)
    int secondTreeSuperArcId = getVertexSuperArcId(vertexId[i], &secondTree);
    // Root Super Arc Id (second tree) of the sub tree containing the vertex and
    // rooted at the value of the vertex in the first scalar field
    int rootSuperArcId = -1;
    if(!isSuperArcAlreadyVisited[secondTreeSuperArcId]) {
      rootSuperArcId = getSubTreeRootSuperArcId(
        &secondTree, secondTreeSuperArcId, vertexValue[i]);
    } else {
      isMandatory = false;
    }

    // Exploration of the sub tree (if not already eliminated)
    if(isMandatory) {
      subTreeSuperArcId.clear();
      std::queue<int> superArcQueue;
      superArcQueue.push(rootSuperArcId);
      while(isMandatory && (!superArcQueue.empty())) {
        int spaId = superArcQueue.front();
        superArcQueue.pop();
        if(isSuperArcAlreadyVisited[spaId]) {
          isMandatory = false;
        } else {
          int downNodeId = secondTree.getSuperArc(spaId)->getDownNodeId();
          int numberOfDownSuperArcs
            = secondTree.getNode(downNodeId)->getNumberOfDownSuperArcs();
          if(numberOfDownSuperArcs > 0) {
            for(int j = 0; j < numberOfDownSuperArcs; j++) {
              superArcQueue.push(
                secondTree.getNode(downNodeId)->getDownSuperArcId(j));
            }
          }
          subTreeSuperArcId.push_back(spaId);
        }
      }
    }

    // If it is a mandatory extremum
    if(isMandatory) {
      // Add it to the list
      mandatoryExtremum.push_back(vertexId[i]);
      // Mark all the super arcs in the sub tree as visited and compute the
      // critical interval
      criticalInterval.push_back(
        std::pair<double, double>(vertexValue[i], vertexValue[i]));
      for(int j = 0; j < (int)subTreeSuperArcId.size(); j++) {
        isSuperArcAlreadyVisited[subTreeSuperArcId[j]] = true;
        int downNodeId
          = secondTree.getSuperArc(subTreeSuperArcId[j])->getDownNodeId();
        double downNodeValue = secondTree.getNodeScalar(downNodeId);
        if(pointType == PointType::Minimum) {
          if(downNodeValue < criticalInterval.back().first) {
            criticalInterval.back().first = downNodeValue;
          }
        } else {
          if(downNodeValue > criticalInterval.back().second) {
            criticalInterval.back().second = downNodeValue;
          }
        }
      }
      // Mark the super arc from the root of the sub tree to the global root
      int upNodeId = secondTree.getSuperArc(rootSuperArcId)->getUpNodeId();
      int spaId = secondTree.getNode(upNodeId)->getUpSuperArcId(0);
      while((spaId != -1) && !(isSuperArcAlreadyVisited[spaId])) {
        isSuperArcAlreadyVisited[spaId] = true;
        upNodeId = secondTree.getSuperArc(spaId)->getUpNodeId();
        spaId = secondTree.getNode(upNodeId)->getUpSuperArcId(0);
      }
    }

    // // List of sub tree super arc ids (second tree)
    // getSubTreeSuperArcIds(secondTree, rootSuperArcId, subTreeSuperArcId);
    // // Check each super arc in the sub tree in the second tree
    // bool isMandatoryCriticalPoint = true;
    // for(int j=0 ; j<(int)subTreeSuperArcId.size() ; j++){
    //   // If one of the super arc has been visited before, the extremum is not
    //   mandatory if(isSuperArcAlreadyVisited[subTreeSuperArcId[j]]){
    //     isMandatoryCriticalPoint = false;
    //     break;
    //   }
    // }
    // // If it is a mandatory extremum
    // if(isMandatoryCriticalPoint){
    //   // Add it to the list
    //   mandatoryExtremum->push_back(vertexId[i]);
    //   // Mark all the super arcs in the sub tree as visited and compute the
    //   critical interval
    //   criticalInterval->push_back(pair<double,double>(vertexValue[i],vertexValue[i]));
    //   for(int j=0 ; j<(int)subTreeSuperArcId.size() ; j++){
    //     isSuperArcAlreadyVisited[subTreeSuperArcId[j]] = true;
    //     int downNodeId =
    //     secondTree->getSuperArc(subTreeSuperArcId[j])->getDownNodeId();
    //     double downNodeValue = secondTree->getNodeScalar(downNodeId);
    //     if(pointType == PointType::Minimum) {
    //       if(downNodeValue < criticalInterval->back().first) {
    //         criticalInterval->back().first = downNodeValue;
    //       }
    //     } else {
    //       if(downNodeValue > criticalInterval->back().second){
    //         criticalInterval->back().second = downNodeValue;
    //       }
    //     }
    //   }
    // }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
  {
    const std::string pt
      = pointType == PointType::Minimum ? "minima" : "maxima";

    this->printMsg("Computed " + std::to_string(mandatoryExtremum.size())
                     + " mandatory " + pt,
                   1.0, t.getElapsedTime(), this->threadNumber_);
    this->printMsg("List of mandatory " + pt, debug::Priority::DETAIL);

    for(size_t i = 0; i < mandatoryExtremum.size(); i++) {
      std::stringstream msg;
      msg << "  -> " << pt << " (" << std::setw(3) << std::right << i << ") ";
      msg << " \t"
          << "Vertex  " << std::setw(9) << std::left << mandatoryExtremum[i];
      msg << " \tInterval  [ " << std::setw(12) << std::right
          << criticalInterval[i].first << " ; " << std::setprecision(9)
          << std::setw(11) << std::left << criticalInterval[i].second << " ]";
      this->printMsg(msg.str(), debug::Priority::DETAIL);
    }
  }

  return 0;
}

int MandatoryCriticalPoints::enumerateMandatorySaddles(
  const PointType pointType,
  SubLevelSetTree &lowerTree,
  SubLevelSetTree &upperTree,
  const std::vector<int> &mandatoryExtremumVertex,
  std::vector<std::pair<int, int>> &mandatorySaddleVertex,
  std::vector<std::vector<int>> &mandatoryMergedExtrema) {

  // Timer
  Timer t;

  // Object to handle requests for common ancestors
  LowestCommonAncestor lowerLca;
  lowerLca.setDebugLevel(debugLevel_);
  LowestCommonAncestor upperLca;
  upperLca.setDebugLevel(debugLevel_);

  /*
  - Add all mandatory extrema in the node list
  - Find all possible saddles
  - Find ancestor of each node
  - Preprocess requests
  - Compute each pair of extrema
  */

  const unsigned int extremaNumber = mandatoryExtremumVertex.size();

  std::vector<int> upperSaddleList;
  std::vector<int> lowerSaddleList;

  std::vector<int> upperTransverse(upperTree.getNumberOfSuperArcs(), -1);
  std::vector<int> lowerTransverse(lowerTree.getNumberOfSuperArcs(), -1);

  std::vector<std::vector<int>> lowerToUpperLinks;
  std::vector<std::vector<int>> upperToLowerLinks;

  std::vector<std::vector<int>> mergedExtrema;

  // NOTE-julien
  // the loop below is not thread-safe.
  // plus, for all tested examples, it turns out to be even faster in serial.
  //   #pragma omp parallel num_threads(threadNumber_)
  {
// Build the list of all saddles joining the extrema
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)mandatoryExtremumVertex.size(); i++) {
      // Super arc in upper and lower trees
      int upperSuperArcId
        = getVertexSuperArcId(mandatoryExtremumVertex[i], &upperTree);
      int lowerSuperArcId
        = getVertexSuperArcId(mandatoryExtremumVertex[i], &lowerTree);
      bool saddleFound;
      bool rootReached;
      bool multipleSaddleFound;
      // Upper Transverse
      saddleFound = false;
      multipleSaddleFound = false;
      rootReached = false;
      do {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperTreeTransverse)
#endif
        {
          // Check if it's the first time the super arc is visited
          if(++upperTransverse[upperSuperArcId]) {
            // If not, a saddle is found
            saddleFound = true;
            // If two or more previous visits, the saddle is already in the list
            if(upperTransverse[upperSuperArcId] > 1) {
              multipleSaddleFound = true;
            }
          }
        }
        if(!saddleFound) {
          // Get the upper super arc
          int upNodeId = upperTree.getSuperArc(upperSuperArcId)->getUpNodeId();
          upperSuperArcId = upperTree.getNode(upNodeId)->getUpSuperArcId(0);
          if(upperSuperArcId == -1) {
            rootReached = true;
          }
        } else {
          if(!multipleSaddleFound) {
            int saddleId
              = upperTree.getSuperArc(upperSuperArcId)->getDownNodeId();
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
            upperSaddleList.push_back(saddleId);
          }
        }
      } while(!(saddleFound || rootReached));
      // Lower Transverse
      saddleFound = false;
      multipleSaddleFound = false;
      rootReached = false;
      do {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerTreeTransverse)
#endif
        {
          // Check if it's the first time the super arc is visited
          if(++lowerTransverse[lowerSuperArcId]) {
            // If not, a saddle is found
            saddleFound = true;
            // If two or more previous visits, the saddle is already in the list
            if(lowerTransverse[lowerSuperArcId] > 1) {
              multipleSaddleFound = true;
            }
          }
        }
        if(!saddleFound) {
          // Get the upper super arc
          int upNodeId = lowerTree.getSuperArc(lowerSuperArcId)->getUpNodeId();
          lowerSuperArcId = lowerTree.getNode(upNodeId)->getUpSuperArcId(0);
          if(lowerSuperArcId == -1) {
            rootReached = true;
          }
        } else {
          if(!multipleSaddleFound) {
            int saddleId
              = lowerTree.getSuperArc(lowerSuperArcId)->getDownNodeId();
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
            lowerSaddleList.push_back(saddleId);
          }
        }
      } while(!(saddleFound || rootReached));
    }
  }

  // NOTE-julien
  // something is not thread-safe above

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
  {
    // Reinitialize the transverse arrays to -1

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)upperTransverse.size(); i++) {
      upperTransverse[i] = -1;
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)lowerTransverse.size(); i++) {
      lowerTransverse[i] = -1;
    }

// Mark the super arcs with the id of the saddle
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)upperSaddleList.size(); i++) {
      int superArcId
        = upperTree.getNode(upperSaddleList[i])->getUpSuperArcId(0);
      upperTransverse[superArcId] = i;
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)lowerSaddleList.size(); i++) {
      int superArcId
        = lowerTree.getNode(lowerSaddleList[i])->getUpSuperArcId(0);
      lowerTransverse[superArcId] = i;
    }

// Create the nodes in the LCA object
#ifdef TTK_ENABLE_OPENMP
#pragma omp sections
#endif
    {
// Upper lca
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      {
        upperLca.addNodes(mandatoryExtremumVertex.size()
                          + upperSaddleList.size());
      }
// Lower lca
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      {
        lowerLca.addNodes(mandatoryExtremumVertex.size()
                          + lowerSaddleList.size());
      }
    }

// For each extrema, find it's first up saddle
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)mandatoryExtremumVertex.size(); i++) {
      int superArcId;
      int lcaSaddleId;
      // Upper Tree
      superArcId = getVertexSuperArcId(mandatoryExtremumVertex[i], &upperTree);
      while((superArcId != -1) && (upperTransverse[superArcId] == -1)) {
        int nodeId = upperTree.getSuperArc(superArcId)->getUpNodeId();
        superArcId = upperTree.getNode(nodeId)->getUpSuperArcId(0);
      }
      if(superArcId != -1) {
        lcaSaddleId = extremaNumber + upperTransverse[superArcId];
        // Ancestor
        upperLca.getNode(i).setAncestor(lcaSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperLca)
#endif
        upperLca.getNode(lcaSaddleId).addSuccessor(i);
      }

      // Lower Tree
      superArcId = getVertexSuperArcId(mandatoryExtremumVertex[i], &lowerTree);
      while((superArcId != -1) && (lowerTransverse[superArcId] == -1)) {
        int nodeId = lowerTree.getSuperArc(superArcId)->getUpNodeId();
        superArcId = lowerTree.getNode(nodeId)->getUpSuperArcId(0);
      }
      if(superArcId != -1) {
        lcaSaddleId = extremaNumber + lowerTransverse[superArcId];
        // Ancestor
        lowerLca.getNode(i).setAncestor(lcaSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerLca)
#endif
        lowerLca.getNode(lcaSaddleId).addSuccessor(i);
      }
    }

// For each saddle, find it's first up saddle
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)upperSaddleList.size(); i++) {
      int superArcId
        = upperTree.getNode(upperSaddleList[i])->getUpSuperArcId(0);
      do {
        int nodeId = upperTree.getSuperArc(superArcId)->getUpNodeId();
        superArcId = upperTree.getNode(nodeId)->getUpSuperArcId(0);
      } while((superArcId != -1) && (upperTransverse[superArcId] == -1));
      if(superArcId != -1) {
        int lcaSuccessorSaddleId = extremaNumber + i;
        int lcaAncestorSaddleId = extremaNumber + upperTransverse[superArcId];
        // Ancestor
        upperLca.getNode(lcaSuccessorSaddleId).setAncestor(lcaAncestorSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperLca)
#endif
        upperLca.getNode(lcaAncestorSaddleId)
          .addSuccessor(lcaSuccessorSaddleId);
      }
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int i = 0; i < (int)lowerSaddleList.size(); i++) {
      int superArcId
        = lowerTree.getNode(lowerSaddleList[i])->getUpSuperArcId(0);
      do {
        int nodeId = lowerTree.getSuperArc(superArcId)->getUpNodeId();
        superArcId = lowerTree.getNode(nodeId)->getUpSuperArcId(0);
      } while((superArcId != -1) && (lowerTransverse[superArcId] == -1));
      if(superArcId != -1) {
        int lcaSuccessorSaddleId = extremaNumber + i;
        int lcaAncestorSaddleId = extremaNumber + lowerTransverse[superArcId];
        // Ancestor
        lowerLca.getNode(lcaSuccessorSaddleId).setAncestor(lcaAncestorSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerLca)
#endif
        lowerLca.getNode(lcaAncestorSaddleId)
          .addSuccessor(lcaSuccessorSaddleId);
      }
    }

// Preprocess lowest common ancestors requests
#ifdef TTK_ENABLE_OPENMP
#pragma omp sections
#endif
    {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      { upperLca.preprocess(); }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      { lowerLca.preprocess(); }
    }

    // Link lists for each thread
    std::vector<std::vector<int>> localLowerToUpperLinks;
    std::vector<std::vector<int>> localUpperToLowerLinks;
    // Merged extrema list for each thread (lower saddles only)
    std::vector<std::vector<int>> localMergedExtrema;
    // Resize of lists
    localLowerToUpperLinks.resize(lowerSaddleList.size());
    localUpperToLowerLinks.resize(upperSaddleList.size());
    localMergedExtrema.resize(lowerSaddleList.size());
    // Number of iterations (triangular loop without diagonal)
    unsigned int kmax = (extremaNumber * (extremaNumber - 1)) / 2;
// Loop over pairs of extrema
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
    for(int k = 0; k < (int)kmax; k++) {
      unsigned int i = k / extremaNumber;
      unsigned int j = k % extremaNumber;
      if(j <= i) {
        i = extremaNumber - i - 2;
        j = extremaNumber - j - 1;
      }
      int uppLCA = upperLca.query(i, j) - extremaNumber;
      int lowLCA = lowerLca.query(i, j) - extremaNumber;
      localLowerToUpperLinks[lowLCA].push_back(uppLCA);
      localUpperToLowerLinks[uppLCA].push_back(lowLCA);
      localMergedExtrema[lowLCA].push_back(i);
      localMergedExtrema[lowLCA].push_back(j);
    }

// Cleaning of duplicates and fusion of vectors
#ifdef TTK_ENABLE_OPENMP
#pragma omp single
#endif
    {
      lowerToUpperLinks.resize(lowerSaddleList.size());
      upperToLowerLinks.resize(upperSaddleList.size());
      mergedExtrema.resize(lowerSaddleList.size());
    }
    // Lower -> Upper
    for(unsigned int i = 0; i < localLowerToUpperLinks.size(); i++) {
      std::vector<int>::iterator newEnd;
      newEnd = unique(
        localLowerToUpperLinks[i].begin(), localLowerToUpperLinks[i].end());
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerToUpperFusion)
#endif
      {
        lowerToUpperLinks[i].insert(
          lowerToUpperLinks[i].end(),
          make_move_iterator(localLowerToUpperLinks[i].begin()),
          make_move_iterator(newEnd));
      }
      localLowerToUpperLinks[i].clear();
    }
    localLowerToUpperLinks.clear();
    // Upper -> Lower
    for(unsigned int i = 0; i < localUpperToLowerLinks.size(); i++) {
      std::vector<int>::iterator newEnd;
      newEnd = unique(
        localUpperToLowerLinks[i].begin(), localUpperToLowerLinks[i].end());
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperToLowerFusion)
#endif
      {
        upperToLowerLinks[i].insert(
          upperToLowerLinks[i].end(),
          make_move_iterator(localUpperToLowerLinks[i].begin()),
          make_move_iterator(newEnd));
      }
      localUpperToLowerLinks[i].clear();
    }
    localUpperToLowerLinks.clear();
    // Merged extrema
    for(unsigned int i = 0; i < localMergedExtrema.size(); i++) {
      std::vector<int>::iterator newEnd;
      std::sort(localMergedExtrema[i].begin(), localMergedExtrema[i].end());
      newEnd
        = unique(localMergedExtrema[i].begin(), localMergedExtrema[i].end());
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(mergedExtremaFusion)
#endif
      {
        mergedExtrema[i].insert(
          mergedExtrema[i].end(),
          make_move_iterator(localMergedExtrema[i].begin()),
          make_move_iterator(newEnd));
      }
      localMergedExtrema[i].clear();
    }
    localMergedExtrema.clear();
  } // End of omp parallel

  // NOTE-julien:
  // something is not thread-safe above
  // the code below is thread-safe

  // Breadth-first search to identify connected components
  std::vector<bool> isLowerVisited(lowerSaddleList.size(), false);
  std::vector<bool> isUpperVisited(upperSaddleList.size(), false);
  std::vector<std::vector<int>> lowComponent;
  std::vector<std::vector<int>> uppComponent;
  for(unsigned int i = 0; i < lowerSaddleList.size(); i++) {
    if(!isLowerVisited[i]) {
      // New component
      lowComponent.push_back(std::vector<int>());
      uppComponent.push_back(std::vector<int>());
      int componentId = static_cast<int>(lowComponent.size()) - 1;
      lowComponent[componentId].push_back(i);
      isLowerVisited[i] = true;
      // Lists of neighbors
      std::vector<int> nextList;
      std::vector<int> currentList;
      currentList.push_back(i);
      // Pointers
      std::vector<bool> *isVisited = &(isUpperVisited);
      std::vector<std::vector<int>> *component = &(uppComponent);
      std::vector<std::vector<int>> *linkList = &(lowerToUpperLinks);
      while(!currentList.empty()) {
        // For each entry in the list
        for(unsigned int j = 0; j < currentList.size(); j++) {
          // For each neighbors
          for(unsigned int k = 0; k < (*linkList)[currentList[j]].size(); k++) {
            // Get the neighbor id;
            int neighbor = (*linkList)[currentList[j]][k];
            // Check if it's already visited
            if(!(*isVisited)[neighbor]) {
              // If not visited, mark it and add it
              (*isVisited)[neighbor] = true;
              (*component)[componentId].push_back(neighbor);
              nextList.push_back(neighbor);
            }
          }
        }
        // Swap pointers and lists
        isVisited = (isVisited == &(isUpperVisited)) ? &(isLowerVisited)
                                                     : &(isUpperVisited);
        component
          = (component == &(uppComponent)) ? &(lowComponent) : &(uppComponent);
        linkList = (linkList == &(lowerToUpperLinks)) ? &(upperToLowerLinks)
                                                      : &(lowerToUpperLinks);
        currentList.swap(nextList);
        nextList.clear();
      }
    }
  }

  // Find pairs of vertices and list of merged extrema
  int numberOfComponents = static_cast<int>(lowComponent.size());
  mandatorySaddleVertex.resize(numberOfComponents, std::pair<int, int>(-1, -1));
  std::vector<std::vector<int>> mandatoryMergedExtrema_tmp;
  mandatoryMergedExtrema_tmp.resize(numberOfComponents, std::vector<int>());

  for(int i = 0; i < numberOfComponents; i++) {
    int nodeId;
    // Find the vertex with minimum value in f-
    nodeId = lowerSaddleList[lowComponent[i][0]];
    mandatorySaddleVertex[i].first = lowerTree.getNode(nodeId)->getVertexId();
    for(unsigned int j = 0; j < lowComponent[i].size(); j++) {
      // First saddle
      int nId = lowerSaddleList[lowComponent[i][j]];
      double nodeScalar = lowerTree.getNodeScalar(nId);
      double refScalar = 0;
      lowerTree.getVertexScalar(mandatorySaddleVertex[i].first, refScalar);
      if(nodeScalar < refScalar) {
        mandatorySaddleVertex[i].first = lowerTree.getNode(nId)->getVertexId();
      }
    }

    for(unsigned int j = 0; j < lowComponent[i].size(); j++) {
      mandatoryMergedExtrema_tmp[i].insert(
        mandatoryMergedExtrema_tmp[i].end(),
        make_move_iterator(mergedExtrema[lowComponent[i][j]].begin()),
        make_move_iterator(mergedExtrema[lowComponent[i][j]].end()));
      mergedExtrema[lowComponent[i][j]].clear();
    }

    // Find the vertex with maximum value in f+
    nodeId = upperSaddleList[uppComponent[i][0]];
    mandatorySaddleVertex[i].second = upperTree.getNode(nodeId)->getVertexId();
    for(unsigned int j = 0; j < uppComponent[i].size(); j++) {
      // First saddle
      int nId = upperSaddleList[uppComponent[i][j]];
      double nodeScalar = upperTree.getNodeScalar(nId);
      double refScalar = 0;
      upperTree.getVertexScalar(mandatorySaddleVertex[i].second, refScalar);
      if(nodeScalar > refScalar) {
        mandatorySaddleVertex[i].second = upperTree.getNode(nId)->getVertexId();
      }
    }
  }

  // Sort the merged extrema list
  std::vector<std::pair<int, int>> order;
  for(unsigned int i = 0; i < mandatoryMergedExtrema_tmp.size(); i++) {
    std::sort(mandatoryMergedExtrema_tmp[i].begin(),
              mandatoryMergedExtrema_tmp[i].end());
    std::vector<int>::iterator newEnd
      = unique(mandatoryMergedExtrema_tmp[i].begin(),
               mandatoryMergedExtrema_tmp[i].end());
    mandatoryMergedExtrema_tmp[i].resize(
      distance(mandatoryMergedExtrema_tmp[i].begin(), newEnd));
    mandatoryMergedExtrema_tmp[i].shrink_to_fit();
    order.push_back(
      std::pair<int, int>(i, mandatoryMergedExtrema_tmp[i].size()));
  }
  // Sort by number of merged extrema
  mandatorySaddleComparaison cmp;
  std::sort(order.begin(), order.end(), cmp);
  // Reorder for output
  mandatoryMergedExtrema.clear();
  mandatoryMergedExtrema.resize(mandatoryMergedExtrema_tmp.size());
  for(unsigned int i = 0; i < order.size(); i++) {
    mandatoryMergedExtrema_tmp[order[i].first].swap(mandatoryMergedExtrema[i]);
  }

// Getting global max and min
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      if(upperMaximumList_.size()) {
        globalMaximumValue_ = upperVertexScalars_[upperMaximumList_[0]];
        for(size_t i = 0; i < upperMaximumList_.size(); i++) {
          if(upperVertexScalars_[upperMaximumList_[i]] > globalMaximumValue_) {
            globalMaximumValue_ = upperVertexScalars_[upperMaximumList_[i]];
          }
        }
      }
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      if(lowerMinimumList_.size()) {
        globalMinimumValue_ = lowerVertexScalars_[lowerMinimumList_[0]];
        for(size_t i = 0; i < lowerMinimumList_.size(); i++) {
          if(lowerVertexScalars_[lowerMinimumList_[i]] < globalMinimumValue_) {
            globalMinimumValue_ = lowerVertexScalars_[lowerMinimumList_[i]];
          }
        }
      }
    }
  }

  const std::string st
    = pointType == PointType::JoinSaddle ? "join saddles" : "split saddles";

  this->printMsg("Computed " + std::to_string(mandatorySaddleVertex.size())
                   + " mandatory " + st,
                 1.0, t.getElapsedTime(), this->threadNumber_);

  this->printMsg("List of " + st + ":", debug::Priority::DETAIL);
  for(size_t i = 0; i < mandatorySaddleVertex.size(); i++) {
    std::stringstream msg;
    msg << "  -> " << st << " (";
    msg << std::setw(3) << std::right << i << ")";
    msg << "  Vertices  <";
    msg << std::setw(9) << std::right << mandatorySaddleVertex[i].first
        << " ; ";
    msg << std::setw(9) << std::left << mandatorySaddleVertex[i].second << ">";
    msg << "  Interval  [";
    double lowerValue = 0;
    double upperValue = 0;
    lowerTree.getVertexScalar(mandatorySaddleVertex[i].first, lowerValue);
    upperTree.getVertexScalar(mandatorySaddleVertex[i].second, upperValue);
    msg << std::setw(12) << std::right << lowerValue << " ; ";
    msg << std::setw(12) << std::left << upperValue << "]";
    this->printMsg(msg.str(), debug::Priority::DETAIL);
  }

  return 0;
}

int MandatoryCriticalPoints::getSubTreeRootSuperArcId(
  const SubLevelSetTree *tree,
  const int &startingSuperArcId,
  const double &targetValue) const {
  // Starting to the input super arc id
  int superArcId = startingSuperArcId;
  // If the super arc exists
  if(superArcId != -1) {
    // Id of the node at the top of the super arc
    int upNodeId = tree->getSuperArc(superArcId)->getUpNodeId();
    // Value of the vertex associated with the node
    double upNodeValue = tree->getNodeScalar(upNodeId);
    // While condition
    int sign = tree->isJoinTree() ? 1 : -1;
    // Climb up the tree
    while(!((sign * targetValue) < (sign * upNodeValue))) {
      int numberOfUpSuperArcs
        = tree->getNode(upNodeId)->getNumberOfUpSuperArcs();
      if(numberOfUpSuperArcs > 0) {
        superArcId = tree->getNode(upNodeId)->getUpSuperArcId(0);
        upNodeId = tree->getSuperArc(superArcId)->getUpNodeId();
        upNodeValue = tree->getNodeScalar(upNodeId);
      } else {
        break;
      }
    }
    return superArcId;
  }
  return -1;
}

int MandatoryCriticalPoints::findCommonAncestorNodeId(
  const SubLevelSetTree *tree,
  const int &vertexId0,
  const int &vertexId1) const {
  int numberOfSuperArcs = tree->getNumberOfSuperArcs();
  std::vector<bool> isSuperArcVisited(numberOfSuperArcs, false);
  int superArcId0 = getVertexSuperArcId(vertexId0, tree);
  int superArcId1 = getVertexSuperArcId(vertexId1, tree);
  int superArcId = superArcId0;
  do {
    isSuperArcVisited[superArcId] = true;
    superArcId = tree->getNode(tree->getSuperArc(superArcId)->getUpNodeId())
                   ->getUpSuperArcId(0);
  } while(superArcId != -1);
  superArcId = superArcId1;
  while(!isSuperArcVisited[superArcId]) {
    superArcId = tree->getNode(tree->getSuperArc(superArcId)->getUpNodeId())
                   ->getUpSuperArcId(0);
  }
  return tree->getSuperArc(superArcId)->getDownNodeId();
}

void MandatoryCriticalPoints::getSubTreeSuperArcIds(
  const SubLevelSetTree *tree,
  const int &rootSuperArcId,
  std::vector<int> &subTreeSuperArcId) const {
  std::vector<int> listOfSuperArcId;
  std::queue<int> superArcIdsToCompute;
  superArcIdsToCompute.push(rootSuperArcId);
  while(!superArcIdsToCompute.empty()) {
    int superArcId = superArcIdsToCompute.front();
    int downNodeId = tree->getSuperArc(superArcId)->getDownNodeId();
    int numberOfUpSuperArcs
      = tree->getNode(downNodeId)->getNumberOfDownSuperArcs();
    if(numberOfUpSuperArcs > 0) {
      for(int i = 0; i < numberOfUpSuperArcs; i++)
        superArcIdsToCompute.push(
          tree->getNode(downNodeId)->getDownSuperArcId(i));
    }
    listOfSuperArcId.push_back(superArcId);
    superArcIdsToCompute.pop();
  }
  listOfSuperArcId.swap(subTreeSuperArcId);
}

int MandatoryCriticalPoints::simplify(
  const double normalizedThreshold,
  const TreeType treeType,
  const std::vector<std::pair<std::pair<int, int>, double>> &extremaSaddlePair,
  const std::vector<std::vector<int>> &mergedExtrema,
  const int numberOfExtrema,
  std::vector<bool> &extremumSimplified,
  std::vector<bool> &saddleSimplified,
  std::vector<int> &extremumParentSaddle,
  std::vector<int> &saddleParentSaddle) const {

  // Simplification threshold value
  double simplificationThreshold;
  if(normalizedThreshold > 1.0)
    simplificationThreshold = (globalMaximumValue_ - globalMinimumValue_);
  else if(normalizedThreshold < 0.0)
    simplificationThreshold = 0.0;
  else
    simplificationThreshold
      = (globalMaximumValue_ - globalMinimumValue_) * normalizedThreshold;

  int extremaSimplifiedNumber = 0;
  int saddlesSimplifiedNumber = 0;

  // First : simplification of extrema
  extremumSimplified.resize(numberOfExtrema);
  fill(extremumSimplified.begin(), extremumSimplified.end(), false);
  for(size_t i = 0; i < extremaSaddlePair.size(); i++) {
    if(extremaSaddlePair[i].second < simplificationThreshold) {
      extremumSimplified[extremaSaddlePair[i].first.first] = true;
      extremaSimplifiedNumber++;
    } else
      break;
  }

  // Second : simplification of saddles
  saddleSimplified.resize(mergedExtrema.size());
  fill(saddleSimplified.begin(), saddleSimplified.end(), false);
  std::vector<int> nonSimplifiedExtremaNumber(mergedExtrema.size(), 0);
  extremumParentSaddle.resize(numberOfExtrema);
  fill(extremumParentSaddle.begin(), extremumParentSaddle.end(), -1);
  saddleParentSaddle.resize(mergedExtrema.size());
  fill(saddleParentSaddle.begin(), saddleParentSaddle.end(), -1);
  for(size_t i = 0; i < mergedExtrema.size(); i++) {
    saddleParentSaddle[i] = i;
    for(size_t j = 0; j < mergedExtrema[i].size(); j++) {
      if(!extremumSimplified[mergedExtrema[i][j]])
        nonSimplifiedExtremaNumber[i]++;
    }
    if(nonSimplifiedExtremaNumber[i] > 1) {
      // Find one non simplified extremum to test
      int extremum = 0;
      while(extremumSimplified[mergedExtrema[i][extremum]]) {
        extremum++;
      }
      if(!(extremumParentSaddle[mergedExtrema[i][extremum]] == -1)) {
        // Find the root
        int rootSaddleId = extremumParentSaddle[mergedExtrema[i][extremum]];
        int lastSaddleId = -1;
        while(rootSaddleId != lastSaddleId) {
          lastSaddleId = rootSaddleId;
          rootSaddleId = saddleParentSaddle[rootSaddleId];
        }
        // Compare the number of merged extrema
        if(!(nonSimplifiedExtremaNumber[i]
             > nonSimplifiedExtremaNumber[rootSaddleId])) {
          saddleSimplified[i] = true;
          saddlesSimplifiedNumber++;
        }
      }
    } else {
      saddleSimplified[i] = true;
      saddlesSimplifiedNumber++;
    }
    if(!saddleSimplified[i]) {
      for(size_t j = 0; j < mergedExtrema[i].size(); j++) {
        int extremaId = mergedExtrema[i][j];
        if(!extremumSimplified[extremaId]) {
          // If the extrema is not already connected, define the parent
          if(extremumParentSaddle[extremaId] == -1) {
            extremumParentSaddle[extremaId] = i;
          } else {
            // Find the root
            int rootSaddleId = extremumParentSaddle[extremaId];
            int lastSaddleId = -1;
            while(rootSaddleId != lastSaddleId) {
              lastSaddleId = rootSaddleId;
              rootSaddleId = saddleParentSaddle[rootSaddleId];
            }
            // Link the root to the proccessed saddle
            saddleParentSaddle[rootSaddleId] = i;
          }
        }
      }
    }
  }

  const std::string pt = treeType == TreeType::JoinTree ? "minima" : "maxima";

  this->printMsg(
    "#Simplified " + pt + ": " + std::to_string(extremaSimplifiedNumber)
    + " (threshold value = " + std::to_string(simplificationThreshold) + ")");

  if(extremaSimplifiedNumber > 0) {
    this->printMsg("List of simplified " + pt + ":", debug::Priority::DETAIL);
    for(size_t i = 0; i < extremumSimplified.size(); i++) {
      if(extremumSimplified[i]) {
        std::stringstream msg;
        msg << "    (" << i << ")";
        this->printMsg(msg.str(), debug::Priority::DETAIL);
      }
    }
  }

  const std::string st
    = treeType == TreeType::JoinTree ? "join saddles" : "split saddles";

  this->printMsg(
    "#Simplified " + st + ": " + std::to_string(saddlesSimplifiedNumber)
    + " (threshold value = " + std::to_string(simplificationThreshold) + ")");

  if(saddlesSimplifiedNumber > 0) {
    this->printMsg("List of simplified " + st + ":", debug::Priority::DETAIL);
    for(size_t i = 0; i < saddleSimplified.size(); i++) {
      if(saddleSimplified[i]) {
        std::stringstream msg;
        msg << "    (" << i << ")";
        this->printMsg(msg.str(), debug::Priority::DETAIL);
      }
    }
  }

  return 0;
}
