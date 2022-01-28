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
///  Progressive Approach used by default
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa ttkScalarFieldCriticalPoints.cpp %for a usage example.

#pragma once

#include <map>

// base code includes
#include <ProgressiveTopology.h>
#include <Triangulation.h>
#include <UnionFind.h>

namespace ttk {

  class ScalarFieldCriticalPoints : public virtual Debug {

  public:
    ScalarFieldCriticalPoints();

    enum class BACKEND { GENERIC = 0, PROGRESSIVE_TOPOLOGY = 1 };

    /**
     * Execute the package.
     * \param offsets Pointer to order field on vertices
     * \param triangulation Triangulation
     * \return Returns 0 upon success, negative values otherwise.
     *
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p offsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    template <class triangulationType = AbstractTriangulation>
    int execute(const SimplexId *const offsets,
                const triangulationType *triangulation);

    template <class triangulationType = AbstractTriangulation>
    int executeLegacy(const SimplexId *const offsets,
                      const triangulationType *triangulation);

    template <class triangulationType = AbstractTriangulation>
    int executeProgressive(const SimplexId *const offsets,
                           const triangulationType *triangulation);

    template <class triangulationType>
    void checkProgressivityRequirement(const triangulationType *triangulation);

    template <class triangulationType = AbstractTriangulation>
    int getNumberOfLowerUpperComponents(const SimplexId vertexId,
                                        const SimplexId *const offsets,
                                        const triangulationType *triangulation,
                                        ttk::SimplexId &lowerComponentNumber,
                                        ttk::SimplexId &upperComponentNumber,
                                        bool &isLowerOnBoundary,
                                        bool &isUpperOnBoundary) const;

    template <class triangulationType = AbstractTriangulation>
    char getCriticalType(const SimplexId &vertexId,
                         const SimplexId *const offsets,
                         const triangulationType *triangulation) const;

    char getCriticalType(const SimplexId &vertexId,
                         const SimplexId *const offsets,
                         const std::vector<std::pair<SimplexId, SimplexId>>
                           &vertexLinkEdgeList) const;

    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

    inline void
      setOutput(std::vector<std::pair<SimplexId, char>> *criticalPoints) {
      criticalPoints_ = criticalPoints;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *triangulation) {
      // pre-condition functions
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
        triangulation->preconditionVertexStars();
        triangulation->preconditionBoundaryVertices();
      }

      setDomainDimension(triangulation->getDimensionality());
      setVertexNumber(triangulation->getNumberOfVertices());
    }

    inline void setVertexLinkEdgeLists(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeList) {
      vertexLinkEdgeLists_ = edgeList;
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

    void displayStats();

  protected:
    int dimension_{};
    SimplexId vertexNumber_{};
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *vertexLinkEdgeLists_{};
    std::vector<std::pair<SimplexId, char>> *criticalPoints_{};

    bool forceNonManifoldCheck{false};

    // progressive
    BACKEND BackEnd{BACKEND::GENERIC};
    ProgressiveTopology progT_{};
    int StartingResolutionLevel{0};
    int StoppingResolutionLevel{-1};
    bool IsResumable{false};
    double TimeLimit{};
  };
} // namespace ttk

// template functions
template <class triangulationType>
int ttk::ScalarFieldCriticalPoints::execute(
  const SimplexId *const offsets, const triangulationType *triangulation) {

  checkProgressivityRequirement(triangulation);

  switch(BackEnd) {

    case BACKEND::PROGRESSIVE_TOPOLOGY:
      this->executeProgressive(offsets, triangulation);
      break;

    case BACKEND::GENERIC:
      this->executeLegacy(offsets, triangulation);
      break;

    default:
      printErr("No method was selected");
  }

  printMsg(ttk::debug::Separator::L1);
  return 0;
}

template <class triangulationType>
int ttk::ScalarFieldCriticalPoints::executeProgressive(
  const SimplexId *const offsets, const triangulationType *triangulation) {

  progT_.setDebugLevel(debugLevel_);
  progT_.setThreadNumber(threadNumber_);
  progT_.setupTriangulation((ttk::ImplicitTriangulation *)triangulation);
  progT_.setStartingResolutionLevel(StartingResolutionLevel);
  progT_.setStoppingResolutionLevel(StoppingResolutionLevel);
  progT_.setTimeLimit(TimeLimit);
  progT_.setIsResumable(IsResumable);
  progT_.setPreallocateMemory(true);

  progT_.computeProgressiveCP(criticalPoints_, offsets);

  displayStats();
  return 0;
}

template <class triangulationType>
int ttk::ScalarFieldCriticalPoints::executeLegacy(
  const SimplexId *const offsets, const triangulationType *triangulation) {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!dimension_) && ((!triangulation) || (triangulation->isEmpty())))
    return -1;
  if((!vertexNumber_) && ((!triangulation) || (triangulation->isEmpty())))
    return -2;
  if((!vertexLinkEdgeLists_) && (!triangulation))
    return -4;
  if(!criticalPoints_)
    return -5;
#endif

  if(triangulation) {
    vertexNumber_ = triangulation->getNumberOfVertices();
    dimension_ = triangulation->getCellVertexNumber(0) - 1;
  }

  printMsg("Extracting critical points...");

  Timer t;

  std::vector<char> vertexTypes(vertexNumber_);

  if(triangulation) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, offsets, triangulation);
    }
  } else if(vertexLinkEdgeLists_) {
    // legacy implementation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, offsets, (*vertexLinkEdgeLists_)[i]);
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

  return 0;
}

template <class triangulationType>
int ttk::ScalarFieldCriticalPoints::getNumberOfLowerUpperComponents(
  const SimplexId vertexId,
  const SimplexId *const offsets,
  const triangulationType *triangulation,
  ttk::SimplexId &lowerComponentNumber,
  ttk::SimplexId &upperComponentNumber,
  bool &isLowerOnBoundary,
  bool &isUpperOnBoundary) const {

  SimplexId neighborNumber = triangulation->getVertexNeighborNumber(vertexId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = 0;
    triangulation->getVertexNeighbor(vertexId, i, neighborId);

    if(offsets[neighborId] < offsets[vertexId]) {
      lowerNeighbors.push_back(neighborId);
      if(dimension_ == 3) {
        if(triangulation->isVertexOnBoundary(neighborId))
          isLowerOnBoundary = true;
      }
    }

    // upper link
    if(offsets[neighborId] > offsets[vertexId]) {
      upperNeighbors.push_back(neighborId);
      if(dimension_ == 3) {
        if(triangulation->isVertexOnBoundary(neighborId))
          isUpperOnBoundary = true;
      }
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    lowerComponentNumber = 0;
    upperComponentNumber = 1;
    return 0;
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    lowerComponentNumber = 1;
    upperComponentNumber = 0;
    return 0;
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

        bool lower0 = offsets[neighborId0] < offsets[vertexId];

        // connect it to everybody except himself and vertexId
        for(SimplexId k = j + 1; k < cellSize; k++) {

          SimplexId neighborId1 = -1;
          triangulation->getCellVertex(cellId, k, neighborId1);

          if((neighborId1 != neighborId0) && (neighborId1 != vertexId)) {

            bool lower1 = offsets[neighborId1] < offsets[vertexId];

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
                (*seeds)[lowerId0] = UnionFind::makeUnion(
                  (*seeds)[lowerId0], (*seeds)[lowerId1]);
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

  lowerComponentNumber = lowerList.size();
  upperComponentNumber = upperList.size();

  return 0;
}

template <class triangulationType>
char ttk::ScalarFieldCriticalPoints::getCriticalType(
  const SimplexId &vertexId,
  const SimplexId *const offsets,
  const triangulationType *triangulation) const {

  bool isLowerOnBoundary = false, isUpperOnBoundary = false;
  SimplexId lowerComponentNumber, upperComponentNumber;
  getNumberOfLowerUpperComponents(vertexId, offsets, triangulation,
                                  lowerComponentNumber, upperComponentNumber,
                                  isLowerOnBoundary, isUpperOnBoundary);

  if(lowerComponentNumber == 0 && upperComponentNumber == 1) {
    return (char)(CriticalType::Local_minimum);
  } else if(lowerComponentNumber == 1 && upperComponentNumber == 0) {
    return (char)(CriticalType::Local_maximum);
  } else if(lowerComponentNumber == 1 && upperComponentNumber == 1) {

    if((dimension_ == 3) && (triangulation->isVertexOnBoundary(vertexId))) {
      // special case of boundary saddles
      if((isUpperOnBoundary) && (!isLowerOnBoundary))
        return (char)(CriticalType::Saddle1);
      if((!isUpperOnBoundary) && (isLowerOnBoundary))
        return (char)(CriticalType::Saddle2);
    }

    // regular point
    return (char)(CriticalType::Regular);
  } else {
    // saddles
    if(dimension_ == 2) {
      if((lowerComponentNumber == 2 && upperComponentNumber == 1)
         || (lowerComponentNumber == 1 && upperComponentNumber == 2)
         || (lowerComponentNumber == 2 && upperComponentNumber == 2)) {
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
      if(lowerComponentNumber == 2 && upperComponentNumber == 1) {
        return (char)(CriticalType::Saddle1);
      } else if(lowerComponentNumber == 1 && upperComponentNumber == 2) {
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

template <class triangulationType>
void ttk::ScalarFieldCriticalPoints::checkProgressivityRequirement(
  const triangulationType *ttkNotUsed(triangulation)) {
  if(BackEnd == BACKEND::PROGRESSIVE_TOPOLOGY) {
    if(std::is_same<triangulationType, ttk::CompactTriangulation>::value) {

      printWrn("CompactTriangulation detected.");
      printWrn("Defaulting to the generic backend.");

      BackEnd = BACKEND::GENERIC;
    } else if(!std::is_same<triangulationType,
                            ttk::ImplicitTriangulation>::value) {

      printWrn("Explicit triangulation detected.");
      printWrn("Defaulting to the generic backend.");

      BackEnd = BACKEND::GENERIC;
    }
  }
}
