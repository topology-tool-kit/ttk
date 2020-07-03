/// \ingroup base
/// \class ttk::ScalarFieldCriticalPoints
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK processing package for the computation of critical points in PL
/// scalar fields defined on PL manifolds.
///
/// This class computes the list of critical points of the input scalar field
/// and classify them according to their type.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Critical points and curvature for embedded polyhedral surfaces" \n
/// Thomas Banchoff \n
/// American Mathematical Monthly, 1970.
///
/// \sa ttkScalarFieldCriticalPoints.cpp %for a usage example.

#pragma once

#include <map>

// base code includes
#include <Triangulation.h>
#include <UnionFind.h>

namespace ttk {

  class ScalarFieldCriticalPoints : public virtual Debug {

  public:
    ScalarFieldCriticalPoints();

    ~ScalarFieldCriticalPoints();

    /// Execute the package.
    /// \param argment Dummy integer argument.
    /// \return Returns 0 upon success, negative values otherwise.
    template <class dataType, class triangulationType = AbstractTriangulation>
    int execute(const dataType *scalarValues,
                const triangulationType *triangulation);

    template <class dataType, class triangulationType = AbstractTriangulation>
    std::pair<SimplexId, SimplexId> getNumberOfLowerUpperComponents(
      const SimplexId vertexId,
      const dataType *scalarValues,
      const triangulationType *triangulation) const;

    template <class dataType, class triangulationType = AbstractTriangulation>
    char getCriticalType(const dataType *scalarValues,
                         const triangulationType *triangulation,
                         const SimplexId &vertexId) const {

      return getCriticalType<dataType, triangulationType>(
        vertexId, scalarValues, triangulation);
    }

    template <class dataType, class triangulationType = AbstractTriangulation>
    char getCriticalType(const SimplexId &vertexId,
                         const dataType *scalarValues,
                         const triangulationType *triangulation) const;

    template <class dataType>
    char getCriticalType(const SimplexId &vertexId,
                         const dataType *scalarValues,
                         const std::vector<std::pair<SimplexId, SimplexId>>
                           &vertexLinkEdgeList) const;

    template <class dataType>
    static bool isSosHigherThan(const SimplexId &offset0,
                                const dataType &value0,
                                const SimplexId &offset1,
                                const dataType &value1) {

      return ((value0 > value1) || ((value0 == value1) && (offset0 > offset1)));
    }

    template <class dataType>
    static bool isSosLowerThan(const SimplexId &offset0,
                               const dataType &value0,
                               const SimplexId &offset1,
                               const dataType &value1) {

      return ((value0 < value1) || ((value0 == value1) && (offset0 < offset1)));
    }

    int setDomainDimension(const int &dimension) {

      dimension_ = dimension;

      return 0;
    }

    int setOutput(std::vector<std::pair<SimplexId, char>> *criticalPoints) {

      criticalPoints_ = criticalPoints;

      return 0;
    }

    int setupTriangulation(AbstractTriangulation *triangulation) {

      // pre-condition functions
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
        triangulation->preconditionVertexStars();
      }

      setDomainDimension(triangulation->getDimensionality());
      setVertexNumber(triangulation->getNumberOfVertices());

      return 0;
    }

    int setSosOffsets(std::vector<SimplexId> *offsets) {

      sosOffsets_ = offsets;

      return 0;
    }

    int setVertexLinkEdgeLists(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeList) {

      vertexLinkEdgeLists_ = edgeList;

      return 0;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    void setNonManifold(const bool b) {
      forceNonManifoldCheck = b;
    }

  protected:
    int dimension_;
    SimplexId vertexNumber_;
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *vertexLinkEdgeLists_;
    std::vector<std::pair<SimplexId, char>> *criticalPoints_;
    std::vector<SimplexId> *sosOffsets_;
    std::vector<SimplexId> localSosOffSets_;

    bool forceNonManifoldCheck;
  };
} // namespace ttk

// template functions
template <class dataType, class triangulationType>
int ttk::ScalarFieldCriticalPoints::execute(
  const dataType *scalarValues, const triangulationType *triangulation) {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!dimension_) && ((!triangulation) || (triangulation->isEmpty())))
    return -1;
  if((!vertexNumber_) && ((!triangulation) || (triangulation->isEmpty())))
    return -2;
  if(!scalarValues)
    return -3;
  if((!vertexLinkEdgeLists_) && (!triangulation))
    return -4;
  if(!criticalPoints_)
    return -5;
#endif

  if(triangulation) {
    vertexNumber_ = triangulation->getNumberOfVertices();
    dimension_ = triangulation->getCellVertexNumber(0) - 1;
  }

  if(!sosOffsets_) {
    // let's use our own local copy
    sosOffsets_ = &localSosOffSets_;
  }
  if((SimplexId)sosOffsets_->size() != vertexNumber_) {
    Timer preProcess;
    sosOffsets_->resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; i++)
      (*sosOffsets_)[i] = i;

    printMsg("Preprocessed " + std::to_string(vertexNumber_) + " offsets.", 1,
             preProcess.getElapsedTime(), 1);
  }

  printMsg(ttk::debug::Separator::L1);

  printMsg("Extracting critical points...");

  Timer t;

  std::vector<char> vertexTypes(vertexNumber_);

  if(triangulation) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, scalarValues, triangulation);
    }
  } else if(vertexLinkEdgeLists_) {
    // legacy implementation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i]
        = getCriticalType(i, scalarValues, (*vertexLinkEdgeLists_)[i]);
    }
  }

  SimplexId minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
            oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;

  // debug msg
  if(debugLevel_ >= (int)debug::Priority::INFO) {
    if(dimension_ == 3) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case(char)(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case(char)(CriticalType::Saddle1):
            oneSaddleNumber++;
            break;

          case(char)(CriticalType::Saddle2):
            twoSaddleNumber++;
            break;

          case(char)(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case(char)(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    } else if(dimension_ == 2) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case(char)(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case(char)(CriticalType::Saddle1):
            saddleNumber++;
            break;

          case(char)(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case(char)(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    }

    {
      std::vector<std::vector<std::string>> stats;
      stats.push_back({"  #Minima", std::to_string(minimumNumber)});
      if(dimension_ == 3) {
        stats.push_back({"  #1-saddles", std::to_string(oneSaddleNumber)});
        stats.push_back({"  #2-saddles", std::to_string(twoSaddleNumber)});
      }
      if(dimension_ == 2) {
        stats.push_back({"  #Saddles", std::to_string(saddleNumber)});
      }
      stats.push_back({"  #Multi-saddles", std::to_string(monkeySaddleNumber)});
      stats.push_back({"  #Maxima", std::to_string(maximumNumber)});

      printMsg(stats);
    }
  }

  // prepare the output
  criticalPoints_->clear();
  criticalPoints_->reserve(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(vertexTypes[i] != (char)(CriticalType::Regular)) {
      criticalPoints_->emplace_back(i, vertexTypes[i]);
    }
  }

  printMsg("Processed " + std::to_string(vertexNumber_) + " vertices", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 0;
}

template <class dataType, class triangulationType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints::getNumberOfLowerUpperComponents(
    const SimplexId vertexId,
    const dataType *scalarValues,
    const triangulationType *triangulation) const {

  SimplexId neighborNumber = triangulation->getVertexNeighborNumber(vertexId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = 0;
    triangulation->getVertexNeighbor(vertexId, i, neighborId);

    if(isSosLowerThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      lowerNeighbors.push_back(neighborId);
    }

    // upper link
    if(isSosHigherThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      upperNeighbors.push_back(neighborId);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    return std::make_pair(0, 1);
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    return std::make_pair(1, 0);
  }

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  SimplexId vertexStarSize = triangulation->getVertexStarNumber(vertexId);

  for(SimplexId i = 0; i < vertexStarSize; i++) {
    SimplexId cellId = 0;
    triangulation->getVertexStar(vertexId, i, cellId);

    SimplexId cellSize = triangulation->getCellVertexNumber(cellId);
    for(SimplexId j = 0; j < cellSize; j++) {
      SimplexId neighborId0 = -1;
      triangulation->getCellVertex(cellId, j, neighborId0);

      if(neighborId0 != vertexId) {
        // we are on the link

        bool lower0 = isSosLowerThan<dataType>(
          (*sosOffsets_)[neighborId0], scalarValues[neighborId0],
          (*sosOffsets_)[vertexId], scalarValues[vertexId]);

        // connect it to everybody except himself and vertexId
        for(SimplexId k = j + 1; k < cellSize; k++) {

          SimplexId neighborId1 = -1;
          triangulation->getCellVertex(cellId, k, neighborId1);

          if((neighborId1 != neighborId0) && (neighborId1 != vertexId)) {

            bool lower1 = isSosLowerThan<dataType>(
              (*sosOffsets_)[neighborId1], scalarValues[neighborId1],
              (*sosOffsets_)[vertexId], scalarValues[vertexId]);

            std::vector<SimplexId> *neighbors = &lowerNeighbors;
            std::vector<UnionFind *> *seeds = &lowerList;

            if(!lower0) {
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }

            if(lower0 == lower1) {
              // connect their union-find sets!
              SimplexId lowerId0 = -1, lowerId1 = -1;
              for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
                if((*neighbors)[l] == neighborId0) {
                  lowerId0 = l;
                }
                if((*neighbors)[l] == neighborId1) {
                  lowerId1 = l;
                }
              }
              if((lowerId0 != -1) && (lowerId1 != -1)) {
                (*seeds)[lowerId0]
                  = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
                (*seeds)[lowerId1] = (*seeds)[lowerId0];
              }
            }
          }
        }
      }
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId)
               + ": lowerLink-#CC=" + std::to_string(lowerList.size())
               + " upperLink-#CC=" + std::to_string(upperList.size()),
             debug::Priority::VERBOSE);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}

template <class dataType, class triangulationType>
char ttk::ScalarFieldCriticalPoints::getCriticalType(
  const SimplexId &vertexId,
  const dataType *scalarValues,
  const triangulationType *triangulation) const {

  SimplexId downValence, upValence;
  std::tie(downValence, upValence) = getNumberOfLowerUpperComponents<dataType>(
    vertexId, scalarValues, triangulation);

  if(downValence == 0 && upValence == 1) {
    return (char)(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return (char)(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return (char)(CriticalType::Regular);
  } else {
    // saddles
    if(dimension_ == 2) {
      if((downValence == 2 && upValence == 1)
         || (downValence == 1 && upValence == 2)
         || (downValence == 2 && upValence == 2)) {
        // regular saddle
        return (char)(CriticalType::Saddle1);
      } else {
        // monkey saddle, saddle + extremum
        return (char)(CriticalType::Degenerate);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to
        // disambiguate boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if(downValence == 2 && upValence == 1) {
        return (char)(CriticalType::Saddle1);
      } else if(downValence == 1 && upValence == 2) {
        return (char)(CriticalType::Saddle2);
      } else {
        // monkey saddle, saddle + extremum
        return (char)(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return (char)(CriticalType::Regular);
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints::getCriticalType(
  const SimplexId &vertexId,
  const dataType *scalarValues,
  const std::vector<std::pair<SimplexId, SimplexId>> &vertexLink) const {

  std::map<SimplexId, SimplexId> global2LowerLink, global2UpperLink;
  std::map<SimplexId, SimplexId>::iterator neighborIt;

  SimplexId lowerCount = 0, upperCount = 0;

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId = vertexLink[i].first;

    // first vertex
    // lower link search
    if(isSosLowerThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }

    // second vertex
    neighborId = vertexLink[i].second;

    // lower link search
    if(isSosLowerThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan<dataType>(
         (*sosOffsets_)[neighborId], scalarValues[neighborId],
         (*sosOffsets_)[vertexId], scalarValues[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
  }

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId) + " lower link ("
               + std::to_string(lowerCount) + " vertices)",
             debug::Priority::VERBOSE);
    printMsg("Vertex #" + std::to_string(vertexId) + " upper link ("
               + std::to_string(upperCount) + " vertices)",
             debug::Priority::VERBOSE);
  }

  if(!lowerCount) {
    // minimum
    return (char)(CriticalType::Local_minimum);
  }
  if(!upperCount) {
    // maximum
    return (char)(CriticalType::Local_maximum);
  }

  // so far 40% of the computation, that's ok.

  // now enumerate the connected components of the lower and upper links
  // NOTE: a breadth first search might be faster than a UF
  // if so, one would need the one-skeleton data structure, not the edge list
  std::vector<UnionFind> lowerSeeds(lowerCount);
  std::vector<UnionFind> upperSeeds(upperCount);
  std::vector<UnionFind *> lowerList(lowerCount);
  std::vector<UnionFind *> upperList(upperCount);
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = &(lowerSeeds[i]);
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = &(upperSeeds[i]);

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId0 = vertexLink[i].first;
    SimplexId neighborId1 = vertexLink[i].second;

    // process the lower link
    if((isSosLowerThan<dataType>(
         (*sosOffsets_)[neighborId0], scalarValues[neighborId0],
         (*sosOffsets_)[vertexId], scalarValues[vertexId]))
       && (isSosLowerThan<dataType>(
         (*sosOffsets_)[neighborId1], scalarValues[neighborId1],
         (*sosOffsets_)[vertexId], scalarValues[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2LowerLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2LowerLink.find(neighborId1);

      lowerList[n0It->second]
        = makeUnion(lowerList[n0It->second], lowerList[n1It->second]);
      lowerList[n1It->second] = lowerList[n0It->second];
    }

    // process the upper link
    if((isSosHigherThan<dataType>(
         (*sosOffsets_)[neighborId0], scalarValues[neighborId0],
         (*sosOffsets_)[vertexId], scalarValues[vertexId]))
       && (isSosHigherThan<dataType>(
         (*sosOffsets_)[neighborId1], scalarValues[neighborId1],
         (*sosOffsets_)[vertexId], scalarValues[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2UpperLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2UpperLink.find(neighborId1);

      upperList[n0It->second]
        = makeUnion(upperList[n0It->second], upperList[n1It->second]);
      upperList[n1It->second] = upperList[n0It->second];
    }
  }

  // let's remove duplicates
  std::vector<UnionFind *>::iterator it;
  // update the UFs if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId)
               + ": lowerLink-#CC=" + std::to_string(lowerList.size())
               + " upperLink-#CC=" + std::to_string(upperList.size()),
             debug::Priority::VERBOSE);
  }

  if((lowerList.size() == 1) && (upperList.size() == 1))
    // regular point
    return (char)(CriticalType::Regular);
  else {
    // saddles
    if(dimension_ == 2) {
      if((lowerList.size() > 2) || (upperList.size() > 2)) {
        // monkey saddle
        return (char)(CriticalType::Degenerate);
      } else {
        // regular saddle
        return (char)(CriticalType::Saddle1);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to disambiguate
        // boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if((lowerList.size() == 2) && (upperList.size() == 1)) {
        return (char)(CriticalType::Saddle1);
      } else if((lowerList.size() == 1) && (upperList.size() == 2)) {
        return (char)(CriticalType::Saddle2);
      } else {
        // monkey saddle
        return (char)(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return (char)(CriticalType::Regular);
}
