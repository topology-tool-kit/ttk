#include <MandatoryCriticalPoints.h>

using namespace std;
using namespace ttk;

MandatoryCriticalPoints::MandatoryCriticalPoints() {
  inputUpperBoundField_ = NULL;
  inputLowerBoundField_ = NULL;
  outputMandatoryMinimum_ = NULL;
  outputMandatoryJoinSaddle_ = NULL;
  outputMandatorySplitSaddle_ = NULL;
  outputMandatoryMaximum_ = NULL;
  vertexNumber_ = 0;
  triangulation_ = NULL;
  normalizedThreshold_ = 0.0;
  upperJoinTree_.setDebugLevel(debugLevel_);
  lowerJoinTree_.setDebugLevel(debugLevel_);
  upperSplitTree_.setDebugLevel(debugLevel_);
  lowerSplitTree_.setDebugLevel(debugLevel_);
}

MandatoryCriticalPoints::~MandatoryCriticalPoints() {
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
  mdtMaxSplitSaddlePair.clear();
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

int MandatoryCriticalPoints::buildSubTrees() {

  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexNumber_ <= 0)
    return -1;
  if((int)upperVertexScalars_.size() != vertexNumber_)
    return -2;
  if((int)lowerVertexScalars_.size() != vertexNumber_)
    return -3;
  if((int)vertexPositions_.size() != vertexNumber_)
    return -4;
  if(!triangulation_
     || (triangulation_->getNumberOfVertices() != vertexNumber_))
    return -5;
  if((int)vertexSoSoffsets_.size() != vertexNumber_)
    return -6;
#endif

  // upperMaximumList_ and lowerMinimumList_ computation (not sorted by function
  // value)
  lowerMinimumList_.clear();
  upperMaximumList_.clear();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < vertexNumber_; i++) {
    bool isLowerMin = true;
    bool isUpperMax = true;
    SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(i);
    for(SimplexId j = 0; j < neighborNumber; j++) {
      SimplexId neighborId;
      triangulation_->getVertexNeighbor(i, j, neighborId);
      if((lowerVertexScalars_[neighborId] < lowerVertexScalars_[i])
         || ((lowerVertexScalars_[neighborId] == lowerVertexScalars_[i])
             && (vertexSoSoffsets_[neighborId] < vertexSoSoffsets_[i])))
        isLowerMin = false;
      if((upperVertexScalars_[neighborId] > upperVertexScalars_[i])
         || ((upperVertexScalars_[neighborId] == upperVertexScalars_[i])
             && (vertexSoSoffsets_[neighborId] > vertexSoSoffsets_[i])))
        isUpperMax = false;
      if(!isUpperMax && !isLowerMin)
        break;
    }
    if(isLowerMin) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
      { lowerMinimumList_.push_back(i); }
    }
    if(isUpperMax) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
      { upperMaximumList_.push_back(i); }
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections num_threads(threadNumber_)
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {upperJoinTree_.setNumberOfVertices(vertexNumber_);
  upperJoinTree_.setVertexScalars(&upperVertexScalars_);
  upperJoinTree_.setVertexPositions(&vertexPositions_);
  upperJoinTree_.setTriangulation(triangulation_);
  upperJoinTree_.setVertexSoSoffsets(&vertexSoSoffsets_);
  upperJoinTree_.buildExtremumList(upperMinimumList_, true);
  upperJoinTree_.build();
}
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
{
  lowerJoinTree_.setNumberOfVertices(vertexNumber_);
  lowerJoinTree_.setVertexScalars(&lowerVertexScalars_);
  lowerJoinTree_.setVertexPositions(&vertexPositions_);
  lowerJoinTree_.setTriangulation(triangulation_);
  lowerJoinTree_.setVertexSoSoffsets(&vertexSoSoffsets_);
  lowerJoinTree_.setMinimumList(lowerMinimumList_);
  lowerJoinTree_.build();
}
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
{
  upperSplitTree_.setNumberOfVertices(vertexNumber_);
  upperSplitTree_.setVertexScalars(&upperVertexScalars_);
  upperSplitTree_.setVertexPositions(&vertexPositions_);
  upperSplitTree_.setTriangulation(triangulation_);
  upperSplitTree_.setVertexSoSoffsets(&vertexSoSoffsets_);
  upperSplitTree_.setMaximumList(upperMaximumList_);
  upperSplitTree_.build();
}
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
{
  lowerSplitTree_.setNumberOfVertices(vertexNumber_);
  lowerSplitTree_.setVertexScalars(&lowerVertexScalars_);
  lowerSplitTree_.setVertexPositions(&vertexPositions_);
  lowerSplitTree_.setTriangulation(triangulation_);
  lowerSplitTree_.setVertexSoSoffsets(&vertexSoSoffsets_);
  lowerSplitTree_.buildExtremumList(lowerMaximumList_, false);
  lowerSplitTree_.build();
}
}
{
  stringstream msg;
  msg << "[MandatoryCriticalPoints] ";
  msg << "4 SubLevelSetTrees computed in ";
  msg << t.getElapsedTime() << " s. (";
  msg << threadNumber_;
  msg << " thread(s)).";
  msg << endl;
  dMsg(cout, msg.str(), timeMsg);
}
return 0;
}

int MandatoryCriticalPoints::buildMandatoryTree(const TreeType treeType) {

  /* Input Variables */
  // Union Find Structure
  const vector<int> *mdtExtremumParentSaddle = (treeType == TreeType::JoinTree)
                                                 ? &(mdtMinimumParentSaddleId_)
                                                 : &(mdtMaximumParentSaddleId_);
  const vector<int> *mdtSaddleParentSaddle
    = (treeType == TreeType::JoinTree) ? &(mdtJoinSaddleParentSaddleId_)
                                       : &(mdtSplitSaddleParentSaddleId_);
  // Simplification flags
  const vector<bool> *isExtremumSimplified = (treeType == TreeType::JoinTree)
                                               ? &(isMdtMinimumSimplified_)
                                               : &(isMdtMaximumSimplified_);
  const vector<bool> *isSaddleSimplified = (treeType == TreeType::JoinTree)
                                             ? &(isMdtJoinSaddleSimplified_)
                                             : &(isMdtSplitSaddleSimplified_);
  // Intervals
  const vector<pair<double, double>> *extremumInterval
    = (treeType == TreeType::JoinTree) ? &(mandatoryMinimumInterval_)
                                       : &(mandatoryMaximumInterval_);
  // Pairs of saddle vertices for the saddle intervals
  const vector<pair<int, int>> *mandatorySaddleVertices
    = (treeType == TreeType::JoinTree) ? &(mandatoryJoinSaddleVertex_)
                                       : &(mandatorySplitSaddleVertex_);
  // Other Informations
  const int extremaNumber = (treeType == TreeType::JoinTree)
                              ? mandatoryMinimumVertex_.size()
                              : mandatoryMaximumVertex_.size();
  const int saddleNumber = (treeType == TreeType::JoinTree)
                             ? mandatoryJoinSaddleVertex_.size()
                             : mandatorySplitSaddleVertex_.size();
  const PointType extremumType = (treeType == TreeType::JoinTree)
                                   ? PointType::Minimum
                                   : PointType::Maximum;
  const PointType saddleType = (treeType == TreeType::JoinTree)
                                 ? PointType::JoinSaddle
                                 : PointType::SplitSaddle;
  const PointType otherExtremumType = (treeType == TreeType::JoinTree)
                                        ? PointType::Maximum
                                        : PointType::Minimum;
  const double globalOtherExtremumValue = (treeType == TreeType::JoinTree)
                                            ? getGlobalMaximum()
                                            : getGlobalMinimum();

  // TODO Check consistency of the input variables here

  /* Output Variables */
  // Graph object
  Graph *mdtTree;
  // Id of the critical component associated to the graph point
  vector<int> *mdtTreePointComponentId;
  // Critical component type associated to the graph point
  vector<PointType> *mdtTreePointType;
  // Lower and Upper values of the interval
  vector<double> *mdtTreePointLowInterval;
  vector<double> *mdtTreePointUpInterval;
  vector<int> *mdtTreeEdgeSwitchable;
  if(treeType == TreeType::JoinTree) {
    mdtTree = &(mdtJoinTree_);
    mdtTreePointComponentId = &(mdtJoinTreePointComponentId_);
    mdtTreePointType = &(mdtJoinTreePointType_);
    mdtTreePointLowInterval = &(mdtJoinTreePointLowInterval_);
    mdtTreePointUpInterval = &(mdtJoinTreePointUpInterval_);
    mdtTreeEdgeSwitchable = &(mdtJoinTreeEdgeSwitchable_);
  } else {
    mdtTree = &(mdtSplitTree_);
    mdtTreePointComponentId = &(mdtSplitTreePointComponentId_);
    mdtTreePointType = &(mdtSplitTreePointType_);
    mdtTreePointLowInterval = &(mdtSplitTreePointLowInterval_);
    mdtTreePointUpInterval = &(mdtSplitTreePointUpInterval_);
    mdtTreeEdgeSwitchable = &(mdtSplitTreeEdgeSwitchable_);
  }

  /* Preliminaries operations */
  mdtTree->clear();
  mdtTreePointComponentId->clear();
  mdtTreePointComponentId->reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointType->clear();
  mdtTreePointType->reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointLowInterval->clear();
  mdtTreePointLowInterval->reserve(extremaNumber + saddleNumber + 1);
  mdtTreePointUpInterval->clear();
  mdtTreePointUpInterval->reserve(extremaNumber + saddleNumber + 1);
  mdtTreeEdgeSwitchable->clear();

  /* Graph Object Building */
  vector<int> saddleGraphVertex(saddleNumber, -1);
  // Create and connect all extrema to their saddle
  for(int i = 0; i < extremaNumber; i++) {
    // If simplified, do nothing and go to the next extremum
    if((*isExtremumSimplified)[i])
      continue;
    // New point in the graph
    int extremumGraphPoint = mdtTree->addVertex();
    mdtTreePointComponentId->push_back(i);
    mdtTreePointType->push_back(extremumType);

    mdtTreePointLowInterval->push_back((*extremumInterval)[i].first);
    mdtTreePointUpInterval->push_back((*extremumInterval)[i].second);
    // Look for the saddle and connect if there is one
    int parentSaddle = (*mdtExtremumParentSaddle)[i];
    // If no parent saddle, end the loop
    if(parentSaddle == -1)
      break;
    // Create the saddle (if not already)
    if(saddleGraphVertex[parentSaddle] == -1) {
      saddleGraphVertex[parentSaddle] = mdtTree->addVertex();
      mdtTreePointComponentId->push_back(parentSaddle);
      mdtTreePointType->push_back(saddleType);
      int lowerVertex = (*mandatorySaddleVertices)[parentSaddle].first;
      int upperVertex = (*mandatorySaddleVertices)[parentSaddle].second;
      mdtTreePointLowInterval->push_back(lowerVertexScalars_[lowerVertex]);
      mdtTreePointUpInterval->push_back(upperVertexScalars_[upperVertex]);
    }
    // Connect the extrema and saddle
    mdtTree->addEdge(saddleGraphVertex[parentSaddle], extremumGraphPoint);
    // Not switchable (saddle -> extremum)
    mdtTreeEdgeSwitchable->push_back(0);
  }

  // If there is no extremum (=> no saddles), do nothing (empty graph)
  if(mdtTree->getNumberOfVertices() == 0) {
    if(debugLevel_ > infoMsg) {
      stringstream msg;
      msg << "[MandatoryCriticalPoints] Mandatory ";
      if(treeType == TreeType::JoinTree) {
        msg << "join tree";
      } else {
        msg << "split tree";
      }
      msg << "empty.";
      msg << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
    return 0;
  }

  // If there is only one extremum (=> no saddles), connect the extremum to the
  // other global extremum
  if(mdtTree->getNumberOfVertices() == 1) {
    mdtTree->addVertex();
    mdtTreePointComponentId->push_back(-1);
    mdtTreePointType->push_back(otherExtremumType);
    mdtTreePointLowInterval->push_back(globalOtherExtremumValue);
    mdtTreePointUpInterval->push_back(globalOtherExtremumValue);
    mdtTree->addEdge(1, 0);
    mdtTreeEdgeSwitchable->push_back(0);
    if(debugLevel_ > infoMsg) {
      stringstream msg;
      msg << "[MandatoryCriticalPoints] Mandatory ";
      if(treeType == TreeType::JoinTree) {
        msg << "join tree";
      } else {
        msg << "split tree";
      }
      msg << " with only one ";
      if(treeType == TreeType::JoinTree) {
        msg << "minimum";
      } else {
        msg << "maximum";
      }
      msg << ".";
      msg << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  } else {
    // Continue to add the remaining saddles
    // Create and connect all remaining saddles
    for(int i = 0; i < saddleNumber; i++) {

      // If simplified, do nothing and go to the next saddle
      if((*isSaddleSimplified)[i]) {
        continue;
      }
      // Create the graph point if not already
      if(saddleGraphVertex[i] == -1) {
        saddleGraphVertex[i] = mdtTree->addVertex();
        mdtTreePointComponentId->push_back(i);
        mdtTreePointType->push_back(saddleType);
        int lowerVertex = (*mandatorySaddleVertices)[i].first;
        int upperVertex = (*mandatorySaddleVertices)[i].second;
        mdtTreePointLowInterval->push_back(lowerVertexScalars_[lowerVertex]);
        mdtTreePointUpInterval->push_back(upperVertexScalars_[upperVertex]);
      }
      // Look for the saddle above and connect if there is one
      int parentSaddle = (*mdtSaddleParentSaddle)[i];
      // If the parent is different from the saddle itself, create it and
      // connect
      if(parentSaddle != i) {
        // Create the graph point if not already
        if(saddleGraphVertex[parentSaddle] == -1) {
          saddleGraphVertex[parentSaddle] = mdtTree->addVertex();
          mdtTreePointComponentId->push_back(parentSaddle);
          mdtTreePointType->push_back(saddleType);
          int lowerVertex = (*mandatorySaddleVertices)[parentSaddle].first;
          int upperVertex = (*mandatorySaddleVertices)[parentSaddle].second;
          mdtTreePointLowInterval->push_back(lowerVertexScalars_[lowerVertex]);
          mdtTreePointUpInterval->push_back(upperVertexScalars_[upperVertex]);
        }
        // Connect the two saddles
        mdtTree->addEdge(saddleGraphVertex[parentSaddle], saddleGraphVertex[i]);
        // Test if switchable : parentSaddle and i
        mdtTreeEdgeSwitchable->push_back(
          areSaddlesSwitchables(treeType, parentSaddle, i));
      } else { // It is the root, create the global extremum to connect with it
        int globalOtherExtremumGraphPoint = mdtTree->addVertex();
        mdtTreePointComponentId->push_back(-1);
        mdtTreePointType->push_back(otherExtremumType);
        mdtTreePointLowInterval->push_back(globalOtherExtremumValue);
        mdtTreePointUpInterval->push_back(globalOtherExtremumValue);
        mdtTree->addEdge(globalOtherExtremumGraphPoint, saddleGraphVertex[i]);
        mdtTreeEdgeSwitchable->push_back(0);
      }
    }
  }

  /* Debug Messages */
  if(debugLevel_ > infoMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] Building of the ";
    if(treeType == TreeType::JoinTree)
      msg << "mandatory join tree";
    else
      msg << "mandatory split tree";
    msg << " : ";
    msg << mdtTree->getNumberOfVertices() << " points, ";
    msg << mdtTree->getNumberOfEdges() << " edges.";
    msg << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  if(debugLevel_ > advancedInfoMsg) {
    stringstream title;
    title << "[MandatoryCriticalPoints] List of ";
    if(treeType == TreeType::JoinTree)
      title << "mandatory join tree";
    else
      title << "mandatory split tree";
    title << " graph points : " << endl;
    dMsg(cout, title.str(), advancedInfoMsg);
    for(int i = 0; i < mdtTree->getNumberOfVertices(); i++) {
      stringstream msg;
      msg << "  (" << setw(3) << right << i << ") ";
      msg << setw(12) << right;
      switch((*mdtTreePointType)[i]) {
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
      msg << "  id = " << setw(3) << (*mdtTreePointComponentId)[i];
      msg << "  I = [ " << (*mdtTreePointLowInterval)[i];
      msg << " ; " << (*mdtTreePointUpInterval)[i];
      msg << " ]";
      msg << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
  }
  if(debugLevel_ > advancedInfoMsg) {
    stringstream title;
    title << "[MandatoryCriticalPoints] List of ";
    if(treeType == TreeType::JoinTree)
      title << "mandatory join tree";
    else
      title << "mandatory split tree";
    title << " graph edges : " << endl;
    dMsg(cout, title.str(), advancedInfoMsg);
    for(int i = 0; i < mdtTree->getNumberOfEdges(); i++) {
      stringstream msg;
      msg << "  (" << setw(3) << right << i << ") ";
      msg << "Points  ";
      msg << setw(3) << right << mdtTree->getEdge(i)->getVertexIdx().first;
      msg << " -> ";
      msg << setw(3) << left << mdtTree->getEdge(i)->getVertexIdx().second;
      msg << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
  }
  return 0;
}

int MandatoryCriticalPoints::buildPairs(const TreeType treeType) {

  /* Input */
  const vector<pair<int, int>> *saddleList = (treeType == TreeType::JoinTree)
                                               ? &(mandatoryJoinSaddleVertex_)
                                               : &(mandatorySplitSaddleVertex_);
  const vector<vector<int>> *mergedExtrema = (treeType == TreeType::JoinTree)
                                               ? &(mergedMinimaId_)
                                               : &(mergedMaximaId_);
  const vector<pair<double, double>> *extremumInterval
    = (treeType == TreeType::JoinTree) ? &(mandatoryMinimumInterval_)
                                       : &(mandatoryMaximumInterval_);
  SubLevelSetTree *lowerTree
    = (treeType == TreeType::JoinTree) ? &(lowerJoinTree_) : &(lowerSplitTree_);
  SubLevelSetTree *upperTree
    = (treeType == TreeType::JoinTree) ? &(upperJoinTree_) : &(upperSplitTree_);
  /* Output */
  vector<pair<pair<int, int>, double>> *extremaSaddlePair
    = (treeType == TreeType::JoinTree) ? &(mdtMinJoinSaddlePair_)
                                       : &(mdtMaxSplitSaddlePair);

#ifndef TTK_ENABLE_KAMIKAZE
  if(lowerTree->isJoinTree() != upperTree->isJoinTree())
    return -1;
#endif

  // List of pairs
  // .first.first = saddle id
  // .first.second = extremum id
  // .second = metric d(M,S)
  extremaSaddlePair->clear();
  // Build list of pairs
  for(int i = 0; i < (int)mergedExtrema->size(); i++) {
    for(int j = 0; j < (int)(*mergedExtrema)[i].size(); j++) {
      // Build a pair (Si,Mj)
      pair<pair<int, int>, double> constructedPair;
      constructedPair.first.first = (*mergedExtrema)[i][j];
      constructedPair.first.second = i;
      // Get the values to compute d(Mj,Si)
      double lowerValue = 0;
      double upperValue = 0;
      if(lowerTree->isJoinTree()) {
        lowerValue = (*extremumInterval)[constructedPair.first.first].first;
        // lowerTree->getVertexScalar(extremumList[(*mergedExtrema)[i][j]],
        // lowerValue);
        upperTree->getVertexScalar((*saddleList)[i].second, upperValue);
      } else {
        lowerTree->getVertexScalar((*saddleList)[i].first, lowerValue);
        // upperTree->getVertexScalar(extremumList[(*mergedExtrema)[i][j]],
        // upperValue);
        upperValue = (*extremumInterval)[constructedPair.first.first].second;
      }
      // Evaluate d(Si,Mj)
      constructedPair.second = fabs(upperValue - lowerValue);
      // Add the pair to the list
      extremaSaddlePair->push_back(constructedPair);
    }
  }
  // Sort pair by increasing value of d(S,M) (.second)
  criticalPointPairComparaison pairComparaison;
  sort(extremaSaddlePair->begin(), extremaSaddlePair->end(), pairComparaison);

  if(debugLevel_ > advancedInfoMsg) {
    stringstream title;
    title << "[MandatoryCriticalPoints] List of ";
    if(lowerTree->isJoinTree())
      title << "minimum - join saddle";
    else
      title << "maximum - split saddle";
    title << " pairs : " << endl;
    dMsg(cout, title.str(), advancedInfoMsg);
    for(int i = 0; i < (int)extremaSaddlePair->size(); i++) {
      stringstream msg;
      msg << "  (" << setw(3) << (*extremaSaddlePair)[i].first.first << ";";
      msg << setw(3) << (*extremaSaddlePair)[i].first.second << ")";
      msg << " -> d = " << (*extremaSaddlePair)[i].second << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
  }
  return 0;
}

int MandatoryCriticalPoints::computePlanarLayout(const TreeType &treeType) {

  /* Input */
  // Mandatory Tree
  const Graph *mdtTree
    = (treeType == TreeType::JoinTree) ? &(mdtJoinTree_) : &(mdtSplitTree_);
  const vector<PointType> *mdtTreePointType = (treeType == TreeType::JoinTree)
                                                ? &(mdtJoinTreePointType_)
                                                : &(mdtSplitTreePointType_);
  const vector<double> *mdtTreePointLowInterval
    = (treeType == TreeType::JoinTree) ? &(mdtJoinTreePointLowInterval_)
                                       : &(mdtSplitTreePointLowInterval_);
  const vector<double> *mdtTreePointUpInterval
    = (treeType == TreeType::JoinTree) ? &(mdtJoinTreePointUpInterval_)
                                       : &(mdtSplitTreePointUpInterval_);

  /* Output */
  // X coordinates
  vector<double> *xCoord = (treeType == TreeType::JoinTree)
                             ? &(mdtJoinTreePointXCoord_)
                             : &(mdtSplitTreePointXCoord_);
  vector<double> *yCoord = (treeType == TreeType::JoinTree)
                             ? &(mdtJoinTreePointYCoord_)
                             : &(mdtSplitTreePointYCoord_);

  /* Informations */
  const int numberOfPoints = mdtTree->getNumberOfVertices();
  const double rangeMin = getGlobalMinimum();
  const double rangeMax = getGlobalMaximum();
  const double range = rangeMax - rangeMin;

  // Get the root
  int rootGraphPointId = -1;
  if(treeType == TreeType::JoinTree) {
    for(size_t i = 0; i < mdtTreePointType->size(); i++) {
      if((*mdtTreePointType)[i] == PointType::Maximum) {
        rootGraphPointId = i;
        break;
      }
    }
  } else {
    for(size_t i = 0; i < mdtTreePointType->size(); i++) {
      if((*mdtTreePointType)[i] == PointType::Minimum) {
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
  if(mdtTree->getVertex(rootGraphPointId)->getNumberOfEdges() == 1) {
    int edgeId = mdtTree->getVertex(rootGraphPointId)->getEdgeIdx(0);
    if(mdtTree->getEdge(edgeId)->getVertexIdx().first == rootGraphPointId) {
      downRootPointId = mdtTree->getEdge(edgeId)->getVertexIdx().second;
    }
  }
  if(downRootPointId == -1) {
    return -2;
  }

  // Resize of outputs
  xCoord->resize(numberOfPoints);
  yCoord->resize(numberOfPoints);

  // Graph Point x intervals
  vector<pair<double, double>> xInterval(numberOfPoints);
  xInterval[rootGraphPointId].first = -0.5 * (double)numberOfPoints;
  xInterval[rootGraphPointId].second = 0.5 * (double)numberOfPoints;
  // Order
  vector<pair<int, double>> xOrder(numberOfPoints);

  queue<int> graphPointQueue;
  graphPointQueue.push(rootGraphPointId);

  while(!graphPointQueue.empty()) {

    int graphPoint = graphPointQueue.front();
    graphPointQueue.pop();

    // Number Of Down Edges
    int numberOfEdges = mdtTree->getVertex(graphPoint)->getNumberOfEdges();
    int numberOfDownEdges = 0;
    for(int i = 0; i < numberOfEdges; i++) {
      int edgeId = mdtTree->getVertex(graphPoint)->getEdgeIdx(i);
      if(mdtTree->getEdge(edgeId)->getVertexIdx().first == graphPoint) {
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
      int edgeId = mdtTree->getVertex(graphPoint)->getEdgeIdx(i);
      int downGraphPoint = -1;
      if(mdtTree->getEdge(edgeId)->getVertexIdx().first == graphPoint) {
        downGraphPoint = mdtTree->getEdge(edgeId)->getVertexIdx().second;
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
    (*yCoord)[i]
      = ((0.5 * ((*mdtTreePointLowInterval)[i] + (*mdtTreePointUpInterval)[i]))
         - rangeMin)
        / range;
  }

  /* X coordinates */
  // Sorting
  pairComparaison xCoordCmp;
  sort(xOrder.begin(), xOrder.end(), xCoordCmp);
  // X coordinates for all points except root (global extremum)
  int pointCount = 0;
  for(int i = 0; i < numberOfPoints; i++) {
    if(xOrder[i].first != rootGraphPointId) {
      (*xCoord)[xOrder[i].first]
        = (double)pointCount * (1.0 / ((int)numberOfPoints - 2.0));
      pointCount++;
    }
  }
  // X coordinate for the root
  (*xCoord)[rootGraphPointId] = (*xCoord)[downRootPointId];

  if(debugLevel_ > advancedInfoMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] Planar layout for mandatory ";
    if(treeType == TreeType::JoinTree)
      msg << "join tree";
    else
      msg << "split tree";
    msg << " : root point " << rootGraphPointId;
    msg << " ( -> " << downRootPointId << ")";
    msg << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
    for(int i = 0; i < numberOfPoints; i++) {
      stringstream msg2;
      msg2 << "  (" << i << ") :"
           << "  x = " << (*xCoord)[i] << "  y = " << (*yCoord)[i];
      msg2 << endl;
      dMsg(cout, msg2.str(), advancedInfoMsg);
    }
  }

  return 0;
}

int MandatoryCriticalPoints::computeExtremumComponent(
  const int &componentId, const PointType &pointType) {

  const SubLevelSetTree *tree = (pointType == PointType::Minimum)
                                  ? &(lowerJoinTree_)
                                  : &(upperSplitTree_);
  const int seedVertexId = (pointType == PointType::Minimum)
                             ? mandatoryMinimumVertex_[componentId]
                             : mandatoryMaximumVertex_[componentId];
  const double value = (pointType == PointType::Minimum)
                         ? upperVertexScalars_[seedVertexId]
                         : lowerVertexScalars_[seedVertexId];
  vector<int> *componentVertexList
    = (pointType == PointType::Minimum)
        ? &(mandatoryMinimumComponentVertices_[componentId])
        : &(mandatoryMaximumComponentVertices_[componentId]);

  // Clear the list
  componentVertexList->clear();
  // Get the super arc id of the vertex
  int superArcId = getVertexSuperArcId(seedVertexId, tree);
  // Get the sub tree root super arc
  int rootSuperArcId = getSubTreeRootSuperArcId(tree, superArcId, value);
  // Get the list of the sub tree super arc ids
  vector<int> subTreeSuperArcId;
  getSubTreeSuperArcIds(tree, rootSuperArcId, subTreeSuperArcId);
  // Comparaison
  int sign = (pointType == PointType::Minimum) ? 1 : -1;
  // Compute each super arc
  for(int i = 0; i < (int)subTreeSuperArcId.size(); i++) {
    const SuperArc *superArc = tree->getSuperArc(subTreeSuperArcId[i]);
    const int numberOfRegularNodes = superArc->getNumberOfRegularNodes();
    // Root super arc and others treated differently
    if(subTreeSuperArcId[i] == rootSuperArcId) {
      // Test the value for each regular node
      for(int j = 0; j < numberOfRegularNodes; j++) {
        int regularNodeId = superArc->getRegularNodeId(j);
        double nodeScalar = tree->getNodeScalar(regularNodeId);
        if(!((sign * nodeScalar) > (sign * value))) {
          int vertexId = tree->getNode(regularNodeId)->getVertexId();
          componentVertexList->push_back(vertexId);
        }
      }
      // Down node
      int downNodeId = superArc->getDownNodeId();
      double nodeScalar = tree->getNodeScalar(downNodeId);
      if(!((sign * nodeScalar) > (sign * value))) {
        int vertexId = tree->getNode(downNodeId)->getVertexId();
        componentVertexList->push_back(vertexId);
      }
    } else {
      // Take all regular nodes
      for(int j = 0; j < numberOfRegularNodes; j++) {
        int regularNodeId = superArc->getRegularNodeId(j);
        int vertexId = tree->getNode(regularNodeId)->getVertexId();
        componentVertexList->push_back(vertexId);
      }
      // Take down node
      int downNodeId = superArc->getDownNodeId();
      int vertexId = tree->getNode(downNodeId)->getVertexId();
      componentVertexList->push_back(vertexId);
    }
  }
  return 0;
}

int MandatoryCriticalPoints::computeSaddleComponent(
  const int &componentId, const PointType &pointType) {

  const int &seedVertexId = (pointType == PointType::JoinSaddle)
                              ? mandatoryJoinSaddleVertex_[componentId].first
                              : mandatorySplitSaddleVertex_[componentId].second;
  const double lowInterval
    = (pointType == PointType::JoinSaddle)
        ? lowerVertexScalars_[mandatoryJoinSaddleVertex_[componentId].first]
        : lowerVertexScalars_[mandatorySplitSaddleVertex_[componentId].first];
  const double upInterval
    = (pointType == PointType::JoinSaddle)
        ? upperVertexScalars_[mandatoryJoinSaddleVertex_[componentId].second]
        : upperVertexScalars_[mandatorySplitSaddleVertex_[componentId].second];
  vector<int> *componentVertexList
    = (pointType == PointType::JoinSaddle)
        ? &(mandatoryJoinSaddleComponentVertices_[componentId])
        : &(mandatorySplitSaddleComponentVertices_[componentId]);

  componentVertexList->clear();

  vector<bool> isVisited(vertexNumber_, false);
  queue<SimplexId> idQueue;
  idQueue.push(seedVertexId);

  while(!(idQueue.empty())) {
    int vertexId = idQueue.front();
    idQueue.pop();
    if(!isVisited[vertexId]) {
      isVisited[vertexId] = true;
      double lowerValue = lowerVertexScalars_[vertexId];
      double upperValue = upperVertexScalars_[vertexId];
      if((pointType == PointType::JoinSaddle && (!(lowerValue > upInterval))
          && (upperValue > lowInterval))
         || (pointType == PointType::SplitSaddle
             && (!(upperValue < lowInterval)) && (lowerValue < upInterval))) {
        componentVertexList->push_back(vertexId);
        // Neighbors
        SimplexId neighborNumber
          = triangulation_->getVertexNeighborNumber(vertexId);
        for(SimplexId i = 0; i < neighborNumber; i++) {
          SimplexId neighborVertexId;
          triangulation_->getVertexNeighbor(vertexId, i, neighborVertexId);
          idQueue.push(neighborVertexId);
        }
      }
    }
  }
  return 0;
}

int MandatoryCriticalPoints::enumerateMandatoryExtrema(
  const PointType pointType) {

  Timer t;

  /* Input */
  SubLevelSetTree *firstTree = (pointType == PointType::Minimum)
                                 ? &(upperJoinTree_)
                                 : &(lowerSplitTree_);
  SubLevelSetTree *secondTree = (pointType == PointType::Minimum)
                                  ? &(lowerJoinTree_)
                                  : &(upperSplitTree_);
  /* Output */
  vector<int> *mandatoryExtremum = (pointType == PointType::Minimum)
                                     ? &(mandatoryMinimumVertex_)
                                     : &(mandatoryMaximumVertex_);
  vector<pair<double, double>> *criticalInterval
    = (pointType == PointType::Minimum) ? &(mandatoryMinimumInterval_)
                                        : &(mandatoryMaximumInterval_);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!(firstTree->isJoinTree() != firstTree->isSplitTree()))
    return -1;
  if(!(secondTree->isJoinTree() != secondTree->isSplitTree()))
    return -2;
  if(!((firstTree->isJoinTree() && secondTree->isJoinTree())
       || (firstTree->isSplitTree() && secondTree->isSplitTree())))
    return -3;
#endif

  // Extremum list in the first tree
  const vector<int> *extremumList = firstTree->getExtremumList();
  // Some tmp variables
  int extremumNumber = (int)extremumList->size();
  vector<int> vertexId(extremumNumber);
  vector<double> vertexValue(extremumNumber);
  vector<int> superArcId(extremumNumber);
  // vector<int> rootSuperArcId(extremumNumber);
  vector<int> subTreeSuperArcId; // vector<vector<int> >
                                 // subTreeSuperArcId(extremumNumber);
  // Clear the list of mandatory extrema
  mandatoryExtremum->clear();
  criticalInterval->clear();
  // To mark the super arc in the second tree as already used for a mandatory
  // critical component
  vector<bool> isSuperArcAlreadyVisited(
    secondTree->getNumberOfSuperArcs(), false);

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  for(int i = 0; i < extremumNumber; i++) {
    // Mandatory until proven otherwise
    bool isMandatory = true;
    // Vertex Id (first tree)
    vertexId[i] = (*extremumList)[i];
    // Vertex Value (first tree)
    firstTree->getVertexScalar(vertexId[i], vertexValue[i]);
    // Super Arc Id (second tree)
    int secondTreeSuperArcId = getVertexSuperArcId(vertexId[i], secondTree);
    // Root Super Arc Id (second tree) of the sub tree containing the vertex and
    // rooted at the value of the vertex in the first scalar field
    int rootSuperArcId = -1;
    if(!isSuperArcAlreadyVisited[secondTreeSuperArcId]) {
      rootSuperArcId = getSubTreeRootSuperArcId(
        secondTree, secondTreeSuperArcId, vertexValue[i]);
    } else {
      isMandatory = false;
    }

    // Exploration of the sub tree (if not already eliminated)
    if(isMandatory) {
      subTreeSuperArcId.clear();
      queue<int> superArcQueue;
      superArcQueue.push(rootSuperArcId);
      while(isMandatory && (!superArcQueue.empty())) {
        int spaId = superArcQueue.front();
        superArcQueue.pop();
        if(isSuperArcAlreadyVisited[spaId]) {
          isMandatory = false;
        } else {
          int downNodeId = secondTree->getSuperArc(spaId)->getDownNodeId();
          int numberOfDownSuperArcs
            = secondTree->getNode(downNodeId)->getNumberOfDownSuperArcs();
          if(numberOfDownSuperArcs > 0) {
            for(int j = 0; j < numberOfDownSuperArcs; j++) {
              superArcQueue.push(
                secondTree->getNode(downNodeId)->getDownSuperArcId(j));
            }
          }
          subTreeSuperArcId.push_back(spaId);
        }
      }
    }

    // If it is a mandatory extremum
    if(isMandatory) {
      // Add it to the list
      mandatoryExtremum->push_back(vertexId[i]);
      // Mark all the super arcs in the sub tree as visited and compute the
      // critical interval
      criticalInterval->push_back(
        pair<double, double>(vertexValue[i], vertexValue[i]));
      for(int j = 0; j < (int)subTreeSuperArcId.size(); j++) {
        isSuperArcAlreadyVisited[subTreeSuperArcId[j]] = true;
        int downNodeId
          = secondTree->getSuperArc(subTreeSuperArcId[j])->getDownNodeId();
        double downNodeValue = secondTree->getNodeScalar(downNodeId);
        if(pointType == PointType::Minimum) {
          if(downNodeValue < criticalInterval->back().first) {
            criticalInterval->back().first = downNodeValue;
          }
        } else {
          if(downNodeValue > criticalInterval->back().second) {
            criticalInterval->back().second = downNodeValue;
          }
        }
      }
      // Mark the super arc from the root of the sub tree to the global root
      int upNodeId = secondTree->getSuperArc(rootSuperArcId)->getUpNodeId();
      int spaId = secondTree->getNode(upNodeId)->getUpSuperArcId(0);
      while((spaId != -1) && !(isSuperArcAlreadyVisited[spaId])) {
        isSuperArcAlreadyVisited[spaId] = true;
        upNodeId = secondTree->getSuperArc(spaId)->getUpNodeId();
        spaId = secondTree->getNode(upNodeId)->getUpSuperArcId(0);
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

  if(debugLevel_ > timeMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] ";
    msg << mandatoryExtremum->size() << " mandatory ";
    if(pointType == PointType::Minimum)
      msg << "minima";
    else
      msg << "maxima";
    msg << " computed in ";
    msg << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))";
    msg << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  if(debugLevel_ > advancedInfoMsg) {
    stringstream title;
    stringstream CPType;
    title << "[MandatoryCriticalPoints] List of mandatory ";
    if(pointType == PointType::Minimum) {
      CPType << "Minimum";
      title << "minima : ";
    } else {
      CPType << "Maximum";
      title << "maxima : ";
    }
    title << endl;
    dMsg(cout, title.str(), advancedInfoMsg);
    for(int i = 0; i < (int)mandatoryExtremum->size(); i++) {
      stringstream msg;
      msg << "  -> " << CPType.str() << " (" << setw(3) << right << i << ") ";
      msg << " \t"
          << "Vertex  " << setw(9) << left << (*mandatoryExtremum)[i];
      msg << " \tInterval  [ " << setw(12) << right
          << (*criticalInterval)[i].first << " ; " << setprecision(9)
          << setw(11) << left << (*criticalInterval)[i].second << " ]";
      msg << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
  }
  return 0;
}

int MandatoryCriticalPoints::enumerateMandatorySaddles(
  const PointType pointType) {
  // Timer
  Timer t;

  /* Input */
  SubLevelSetTree *lowerTree = (pointType == PointType::JoinSaddle)
                                 ? &(lowerJoinTree_)
                                 : &(lowerSplitTree_);
  SubLevelSetTree *upperTree = (pointType == PointType::JoinSaddle)
                                 ? &(upperJoinTree_)
                                 : &(upperSplitTree_);
  const vector<int> *mandatoryExtremumVertex
    = (pointType == PointType::JoinSaddle) ? &(mandatoryMinimumVertex_)
                                           : &(mandatoryMaximumVertex_);
  /* Output */
  vector<pair<int, int>> *mandatorySaddleVertex
    = (pointType == PointType::JoinSaddle) ? &(mandatoryJoinSaddleVertex_)
                                           : &(mandatorySplitSaddleVertex_);
  vector<vector<int>> *mandatoryMergedExtrema
    = (pointType == PointType::JoinSaddle) ? &(mergedMinimaId_)
                                           : &(mergedMaximaId_);

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

  unsigned int extremaNumber = mandatoryExtremumVertex->size();

  vector<int> upperSaddleList;
  vector<int> lowerSaddleList;

  vector<int> upperTransverse(upperTree->getNumberOfSuperArcs(), -1);
  vector<int> lowerTransverse(lowerTree->getNumberOfSuperArcs(), -1);

  vector<vector<int>> lowerToUpperLinks;
  vector<vector<int>> upperToLowerLinks;

  vector<vector<int>> mergedExtrema;

  // NOTE-julien
  // the loop below is not thread-safe.
  // plus, for all tested examples, it turns out to be even faster in serial.
  //   #pragma omp parallel num_threads(threadNumber_)
  {
// Build the list of all saddles joining the extrema
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)mandatoryExtremumVertex->size(); i++) {
      // Super arc in upper and lower trees
      int upperSuperArcId
        = getVertexSuperArcId((*mandatoryExtremumVertex)[i], upperTree);
      int lowerSuperArcId
        = getVertexSuperArcId((*mandatoryExtremumVertex)[i], lowerTree);
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
          int upNodeId = upperTree->getSuperArc(upperSuperArcId)->getUpNodeId();
          upperSuperArcId = upperTree->getNode(upNodeId)->getUpSuperArcId(0);
          if(upperSuperArcId == -1) {
            rootReached = true;
          }
        } else {
          if(!multipleSaddleFound) {
            int saddleId
              = upperTree->getSuperArc(upperSuperArcId)->getDownNodeId();
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
          int upNodeId = lowerTree->getSuperArc(lowerSuperArcId)->getUpNodeId();
          lowerSuperArcId = lowerTree->getNode(upNodeId)->getUpSuperArcId(0);
          if(lowerSuperArcId == -1) {
            rootReached = true;
          }
        } else {
          if(!multipleSaddleFound) {
            int saddleId
              = lowerTree->getSuperArc(lowerSuperArcId)->getDownNodeId();
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
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)upperTransverse.size(); i++) {
      upperTransverse[i] = -1;
    }

#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)lowerTransverse.size(); i++) {
      lowerTransverse[i] = -1;
    }

// Mark the super arcs with the id of the saddle
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)upperSaddleList.size(); i++) {
      int superArcId
        = upperTree->getNode(upperSaddleList[i])->getUpSuperArcId(0);
      upperTransverse[superArcId] = i;
    }

#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)lowerSaddleList.size(); i++) {
      int superArcId
        = lowerTree->getNode(lowerSaddleList[i])->getUpSuperArcId(0);
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
        upperLca.addNodes(mandatoryExtremumVertex->size()
                          + upperSaddleList.size());
      }
// Lower lca
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      {
        lowerLca.addNodes(mandatoryExtremumVertex->size()
                          + lowerSaddleList.size());
      }
    }

// For each extrema, find it's first up saddle
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)mandatoryExtremumVertex->size(); i++) {
      int superArcId;
      int lcaSaddleId;
      // Upper Tree
      superArcId
        = getVertexSuperArcId((*mandatoryExtremumVertex)[i], upperTree);
      while((superArcId != -1) && (upperTransverse[superArcId] == -1)) {
        int nodeId = upperTree->getSuperArc(superArcId)->getUpNodeId();
        superArcId = upperTree->getNode(nodeId)->getUpSuperArcId(0);
      }
      if(superArcId != -1) {
        lcaSaddleId = extremaNumber + upperTransverse[superArcId];
        // Ancestor
        upperLca.getNode(i)->setAncestor(lcaSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperLca)
#endif
        upperLca.getNode(lcaSaddleId)->addSuccessor(i);
      }

      // Lower Tree
      superArcId
        = getVertexSuperArcId((*mandatoryExtremumVertex)[i], lowerTree);
      while((superArcId != -1) && (lowerTransverse[superArcId] == -1)) {
        int nodeId = lowerTree->getSuperArc(superArcId)->getUpNodeId();
        superArcId = lowerTree->getNode(nodeId)->getUpSuperArcId(0);
      }
      if(superArcId != -1) {
        lcaSaddleId = extremaNumber + lowerTransverse[superArcId];
        // Ancestor
        lowerLca.getNode(i)->setAncestor(lcaSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerLca)
#endif
        lowerLca.getNode(lcaSaddleId)->addSuccessor(i);
      }
    }

// For each saddle, find it's first up saddle
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)upperSaddleList.size(); i++) {
      int superArcId
        = upperTree->getNode(upperSaddleList[i])->getUpSuperArcId(0);
      do {
        int nodeId = upperTree->getSuperArc(superArcId)->getUpNodeId();
        superArcId = upperTree->getNode(nodeId)->getUpSuperArcId(0);
      } while((superArcId != -1) && (upperTransverse[superArcId] == -1));
      if(superArcId != -1) {
        int lcaSuccessorSaddleId = extremaNumber + i;
        int lcaAncestorSaddleId = extremaNumber + upperTransverse[superArcId];
        // Ancestor
        upperLca.getNode(lcaSuccessorSaddleId)
          ->setAncestor(lcaAncestorSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(upperLca)
#endif
        upperLca.getNode(lcaAncestorSaddleId)
          ->addSuccessor(lcaSuccessorSaddleId);
      }
    }

#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
#endif
    for(int i = 0; i < (int)lowerSaddleList.size(); i++) {
      int superArcId
        = lowerTree->getNode(lowerSaddleList[i])->getUpSuperArcId(0);
      do {
        int nodeId = lowerTree->getSuperArc(superArcId)->getUpNodeId();
        superArcId = lowerTree->getNode(nodeId)->getUpSuperArcId(0);
      } while((superArcId != -1) && (lowerTransverse[superArcId] == -1));
      if(superArcId != -1) {
        int lcaSuccessorSaddleId = extremaNumber + i;
        int lcaAncestorSaddleId = extremaNumber + lowerTransverse[superArcId];
        // Ancestor
        lowerLca.getNode(lcaSuccessorSaddleId)
          ->setAncestor(lcaAncestorSaddleId);
// Successor
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(lowerLca)
#endif
        lowerLca.getNode(lcaAncestorSaddleId)
          ->addSuccessor(lcaSuccessorSaddleId);
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
    vector<vector<int>> localLowerToUpperLinks;
    vector<vector<int>> localUpperToLowerLinks;
    // Merged extrema list for each thread (lower saddles only)
    vector<vector<int>> localMergedExtrema;
    // Resize of lists
    localLowerToUpperLinks.resize(lowerSaddleList.size());
    localUpperToLowerLinks.resize(upperSaddleList.size());
    localMergedExtrema.resize(lowerSaddleList.size());
    // Number of iterations (triangular loop without diagonal)
    unsigned int kmax = (extremaNumber * (extremaNumber - 1)) / 2;
// Loop over pairs of extrema
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp for
#else
#pragma omp for schedule(auto)
#endif
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
      vector<int>::iterator newEnd;
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
      vector<int>::iterator newEnd;
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
      vector<int>::iterator newEnd;
      sort(localMergedExtrema[i].begin(), localMergedExtrema[i].end());
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
  vector<bool> isLowerVisited(lowerSaddleList.size(), false);
  vector<bool> isUpperVisited(upperSaddleList.size(), false);
  vector<vector<int>> lowComponent;
  vector<vector<int>> uppComponent;
  for(unsigned int i = 0; i < lowerSaddleList.size(); i++) {
    if(!isLowerVisited[i]) {
      // New component
      lowComponent.push_back(vector<int>());
      uppComponent.push_back(vector<int>());
      int componentId = static_cast<int>(lowComponent.size()) - 1;
      lowComponent[componentId].push_back(i);
      isLowerVisited[i] = true;
      // Lists of neighbors
      vector<int> nextList;
      vector<int> currentList;
      currentList.push_back(i);
      // Pointers
      vector<bool> *isVisited = &(isUpperVisited);
      vector<vector<int>> *component = &(uppComponent);
      vector<vector<int>> *linkList = &(lowerToUpperLinks);
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
  (*mandatorySaddleVertex).resize(numberOfComponents, pair<int, int>(-1, -1));
  vector<vector<int>> mandatoryMergedExtrema_tmp;
  mandatoryMergedExtrema_tmp.resize(numberOfComponents, vector<int>());

  for(int i = 0; i < numberOfComponents; i++) {
    int nodeId;
    // Find the vertex with minimum value in f-
    nodeId = lowerSaddleList[lowComponent[i][0]];
    (*mandatorySaddleVertex)[i].first
      = lowerTree->getNode(nodeId)->getVertexId();
    for(unsigned int j = 0; j < lowComponent[i].size(); j++) {
      // First saddle
      int nId = lowerSaddleList[lowComponent[i][j]];
      double nodeScalar = lowerTree->getNodeScalar(nId);
      double refScalar = 0;
      lowerTree->getVertexScalar((*mandatorySaddleVertex)[i].first, refScalar);
      if(nodeScalar < refScalar) {
        (*mandatorySaddleVertex)[i].first
          = lowerTree->getNode(nId)->getVertexId();
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
    (*mandatorySaddleVertex)[i].second
      = upperTree->getNode(nodeId)->getVertexId();
    for(unsigned int j = 0; j < uppComponent[i].size(); j++) {
      // First saddle
      int nId = upperSaddleList[uppComponent[i][j]];
      double nodeScalar = upperTree->getNodeScalar(nId);
      double refScalar = 0;
      upperTree->getVertexScalar((*mandatorySaddleVertex)[i].second, refScalar);
      if(nodeScalar > refScalar) {
        (*mandatorySaddleVertex)[i].second
          = upperTree->getNode(nId)->getVertexId();
      }
    }
  }

  // Sort the merged extrema list
  vector<pair<int, int>> order;
  for(unsigned int i = 0; i < mandatoryMergedExtrema_tmp.size(); i++) {
    sort(mandatoryMergedExtrema_tmp[i].begin(),
         mandatoryMergedExtrema_tmp[i].end());
    vector<int>::iterator newEnd = unique(mandatoryMergedExtrema_tmp[i].begin(),
                                          mandatoryMergedExtrema_tmp[i].end());
    mandatoryMergedExtrema_tmp[i].resize(
      distance(mandatoryMergedExtrema_tmp[i].begin(), newEnd));
    mandatoryMergedExtrema_tmp[i].shrink_to_fit();
    order.push_back(pair<int, int>(i, mandatoryMergedExtrema_tmp[i].size()));
  }
  // Sort by number of merged extrema
  mandatorySaddleComparaison cmp;
  sort(order.begin(), order.end(), cmp);
  // Reorder for output
  (*mandatoryMergedExtrema).clear();
  (*mandatoryMergedExtrema).resize(mandatoryMergedExtrema_tmp.size());
  for(unsigned int i = 0; i < order.size(); i++) {
    mandatoryMergedExtrema_tmp[order[i].first].swap(
      (*mandatoryMergedExtrema)[i]);
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

  // Debug messages
  if(debugLevel_ > timeMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] ";
    msg << mandatorySaddleVertex->size();
    if(pointType == PointType::JoinSaddle)
      msg << " mandatory join saddles computed in ";
    else
      msg << " mandatory split saddles computed in ";
    msg << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))";
    msg << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  if(debugLevel_ > advancedInfoMsg) {
    stringstream title;
    stringstream pointTypeMsg;
    title << "[MandatoryCriticalPoints] List of ";
    if(lowerTree->isJoinTree()) {
      title << "join saddles : ";
      pointTypeMsg << "Join Saddle";
    } else {
      title << "split saddles : ";
      pointTypeMsg << "Split Saddle";
    }
    title << endl;
    dMsg(cout, title.str(), advancedInfoMsg);
    for(int i = 0; i < (int)mandatorySaddleVertex->size(); i++) {
      stringstream msg;
      msg << "  -> " << pointTypeMsg.str() << " (";
      msg << setw(3) << right << i << ")";
      msg << "  Vertices  <";
      msg << setw(9) << right << (*mandatorySaddleVertex)[i].first << " ; ";
      msg << setw(9) << left << (*mandatorySaddleVertex)[i].second << ">";
      msg << "  Interval  [";
      double lowerValue = 0;
      double upperValue = 0;
      lowerTree->getVertexScalar((*mandatorySaddleVertex)[i].first, lowerValue);
      upperTree->getVertexScalar(
        (*mandatorySaddleVertex)[i].second, upperValue);
      msg << setw(12) << right << lowerValue << " ; ";
      msg << setw(12) << left << upperValue << "]";
      msg << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
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
  vector<bool> isSuperArcVisited(numberOfSuperArcs, false);
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
  vector<int> &subTreeSuperArcId) const {
  vector<int> listOfSuperArcId;
  queue<int> superArcIdsToCompute;
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

int MandatoryCriticalPoints::simplify(const double &normalizedThreshold,
                                      const TreeType treeType) {

  /* Input */
  const vector<pair<pair<int, int>, double>> *extremaSaddlePair
    = (treeType == TreeType::JoinTree) ? &(mdtMinJoinSaddlePair_)
                                       : &(mdtMaxSplitSaddlePair);
  const vector<vector<int>> *mergedExtrema = (treeType == TreeType::JoinTree)
                                               ? &(mergedMinimaId_)
                                               : &(mergedMaximaId_);
  const int numberOfExtrema = (treeType == TreeType::JoinTree)
                                ? (int)mandatoryMinimumVertex_.size()
                                : (int)mandatoryMaximumVertex_.size();
  /* Output */
  vector<bool> *extremumSimplified = (treeType == TreeType::JoinTree)
                                       ? &(isMdtMinimumSimplified_)
                                       : &(isMdtMaximumSimplified_);
  vector<bool> *saddleSimplified = (treeType == TreeType::JoinTree)
                                     ? &(isMdtJoinSaddleSimplified_)
                                     : &(isMdtSplitSaddleSimplified_);
  vector<int> *extremumParentSaddle = (treeType == TreeType::JoinTree)
                                        ? &(mdtMinimumParentSaddleId_)
                                        : &(mdtMaximumParentSaddleId_);
  vector<int> *saddleParentSaddle = (treeType == TreeType::JoinTree)
                                      ? &(mdtJoinSaddleParentSaddleId_)
                                      : &(mdtSplitSaddleParentSaddleId_);

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
  extremumSimplified->resize(numberOfExtrema);
  fill(extremumSimplified->begin(), extremumSimplified->end(), false);
  for(int i = 0; i < (int)extremaSaddlePair->size(); i++) {
    if((*extremaSaddlePair)[i].second < simplificationThreshold) {
      (*extremumSimplified)[(*extremaSaddlePair)[i].first.first] = true;
      extremaSimplifiedNumber++;
    } else
      break;
  }

  // Second : simplification of saddles
  saddleSimplified->resize(mergedExtrema->size());
  fill(saddleSimplified->begin(), saddleSimplified->end(), false);
  vector<int> nonSimplifiedExtremaNumber(mergedExtrema->size(), 0);
  extremumParentSaddle->resize(numberOfExtrema);
  fill(extremumParentSaddle->begin(), extremumParentSaddle->end(), -1);
  saddleParentSaddle->resize(mergedExtrema->size());
  fill(saddleParentSaddle->begin(), saddleParentSaddle->end(), -1);
  for(int i = 0; i < (int)mergedExtrema->size(); i++) {
    (*saddleParentSaddle)[i] = i;
    for(int j = 0; j < (int)(*mergedExtrema)[i].size(); j++) {
      if(!(*extremumSimplified)[(*mergedExtrema)[i][j]])
        nonSimplifiedExtremaNumber[i]++;
    }
    if(nonSimplifiedExtremaNumber[i] > 1) {
      // Find one non simplified extremum to test
      int extremum = 0;
      while((*extremumSimplified)[(*mergedExtrema)[i][extremum]]) {
        extremum++;
      }
      if(!((*extremumParentSaddle)[(*mergedExtrema)[i][extremum]] == -1)) {
        // Find the root
        int rootSaddleId
          = (*extremumParentSaddle)[(*mergedExtrema)[i][extremum]];
        int lastSaddleId = -1;
        while(rootSaddleId != lastSaddleId) {
          lastSaddleId = rootSaddleId;
          rootSaddleId = (*saddleParentSaddle)[rootSaddleId];
        }
        // Compare the number of merged extrema
        if(!(nonSimplifiedExtremaNumber[i]
             > nonSimplifiedExtremaNumber[rootSaddleId])) {
          (*saddleSimplified)[i] = true;
          saddlesSimplifiedNumber++;
        }
      }
    } else {
      (*saddleSimplified)[i] = true;
      saddlesSimplifiedNumber++;
    }
    if(!(*saddleSimplified)[i]) {
      for(int j = 0; j < (int)(*mergedExtrema)[i].size(); j++) {
        int extremaId = (*mergedExtrema)[i][j];
        if(!(*extremumSimplified)[extremaId]) {
          // If the extrema is not already connected, define the parent
          if((*extremumParentSaddle)[extremaId] == -1) {
            (*extremumParentSaddle)[extremaId] = i;
          } else {
            // Find the root
            int rootSaddleId = (*extremumParentSaddle)[extremaId];
            int lastSaddleId = -1;
            while(rootSaddleId != lastSaddleId) {
              lastSaddleId = rootSaddleId;
              rootSaddleId = (*saddleParentSaddle)[rootSaddleId];
            }
            // Link the root to the proccessed saddle
            (*saddleParentSaddle)[rootSaddleId] = i;
          }
        }
      }
    }
  }
  /* Debug messages */
  if(debugLevel_ > infoMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] Number of ";
    if(treeType == TreeType::JoinTree)
      msg << "minima";
    else
      msg << "maxima";
    msg << " simplified : " << extremaSimplifiedNumber;
    msg << " (threshold value = " << simplificationThreshold << ")";
    msg << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  if(debugLevel_ > advancedInfoMsg && extremaSimplifiedNumber > 0) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] List of simplified ";
    if(treeType == TreeType::JoinTree)
      msg << "minima";
    else
      msg << "maxima";
    msg << " : " << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
    for(size_t i = 0; i < extremumSimplified->size(); i++) {
      if((*extremumSimplified)[i]) {
        stringstream msg2;
        msg2 << "    (" << i << ")" << endl;
        dMsg(cout, msg2.str(), advancedInfoMsg);
      }
    }
  }
  if(debugLevel_ > infoMsg) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] Number of ";
    if(treeType == TreeType::JoinTree)
      msg << "join saddles";
    else
      msg << "split saddles";
    msg << " simplified : " << saddlesSimplifiedNumber;
    msg << " (threshold value = " << simplificationThreshold << ")";
    msg << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  if(debugLevel_ > advancedInfoMsg && saddlesSimplifiedNumber > 0) {
    stringstream msg;
    msg << "[MandatoryCriticalPoints] List of simplified ";
    if(treeType == TreeType::JoinTree)
      msg << "join saddles";
    else
      msg << "split saddles";
    msg << " : " << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
    for(size_t i = 0; i < saddleSimplified->size(); i++) {
      if((*saddleSimplified)[i]) {
        stringstream msg2;
        msg2 << "    (" << i << ")" << endl;
        dMsg(cout, msg2.str(), advancedInfoMsg);
      }
    }
  }

  return 0;
}
