/// \ingroup base
/// \class ttk::MandatoryCriticalPoints
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2016.
///
/// \brief TTK processing package for the computation of mandatory critical
/// points in uncertain scalar data.
///
/// This package computes the mandatory critical points of uncertain scalar
/// fields defined on PL (nD) manifolds. The input uncertain data is represented
/// by reliable bound fields for each vertex.
///
/// \warning SimplexId (large large datasets). This class builds and runs
/// with the new triangulation API (SimplexId) but may need adjustments when
/// addressing more than integers (large datasets).
///
/// \b Related \b publication \n
/// "Mandatory Critical Points of 2D Uncertain Scalar Fields" \n
/// David Guenther, Joseph Salmon, Julien Tierny \n
/// Proc. of EuroVis 2014. \n
/// Computer Graphics Forum, 2014.
///
/// \sa ttkMandatoryCriticalPoints.cpp %for a usage example.

#pragma once

// base code includes
#include <ContourTree.h>
#include <LowestCommonAncestor.h>
#include <Triangulation.h>

// std includes
#include <queue>
#include <vector>

namespace ttk {

  class Graph : virtual public Debug {
  public:
    class Vertex {
      friend class Graph;

    public:
      inline int getNumberOfEdges() const {
        return (int)edgeIdx_.size();
      }
      inline int getEdgeIdx(const int &connectedEdge) const {
        if(connectedEdge < (int)edgeIdx_.size())
          return edgeIdx_[connectedEdge];
        else
          return -1;
      }

    protected:
      std::vector<int> edgeIdx_{};
    };
    class Edge {
      friend class Graph;

    public:
      Edge(const int &start, const int &end) {
        vertexIdx_.first = start;
        vertexIdx_.second = end;
      }
      const std::pair<int, int> &getVertexIdx() const {
        return vertexIdx_;
      }

    protected:
      std::pair<int, int> vertexIdx_{};
    };

  public:
    inline int getNumberOfVertices() const {
      return (int)vertexList_.size();
    }
    inline int getNumberOfEdges() const {
      return (int)edgeList_.size();
    }
    /// Add a vertex, returns it's index
    int addVertex() {
      Vertex newVertex;
      vertexList_.push_back(newVertex);
      return (int)vertexList_.size() - 1;
    }
    /// Get a pointer to the vertex idx
    const Vertex &getVertex(const int idx) const {
      return vertexList_[idx];
    }
    /// Add an edge between the vertex start and end, returns it's index
    int addEdge(const int &start, const int &end) {
      if((start < (int)vertexList_.size()) && (end < (int)vertexList_.size())) {
        Edge newEdge(start, end);
        edgeList_.push_back(newEdge);
        vertexList_[start].edgeIdx_.push_back(edgeList_.size() - 1);
        vertexList_[end].edgeIdx_.push_back(edgeList_.size() - 1);
        return (int)edgeList_.size() - 1;
      } else
        return -1;
    }
    /// Get a pointer to the edge idx, returns NULL if the idx is incorrect or
    /// if the edge has been removed.
    const Edge &getEdge(const int idx) const {
      return edgeList_[idx];
    }
    inline void clear() {
      vertexList_.clear();
      edgeList_.clear();
    }

  protected:
    std::vector<Vertex> vertexList_{};
    std::vector<Edge> edgeList_{};
  };

  /// Comparison between critical point pairs ( (Extremum,Saddle), dist(M,S) )
  struct criticalPointPairComparaison {
    bool operator()(const std::pair<std::pair<int, int>, double> &left,
                    const std::pair<std::pair<int, int>, double> &right) const {
      return (left.second < right.second);
    }
  };

  /// Comparison between mandatory saddles (Saddle id, Number of merged extrema)
  struct mandatorySaddleComparaison {
    bool operator()(const std::pair<int, int> &left,
                    const std::pair<int, int> &right) const {
      return (left.second < right.second);
    }
  };

  // TODO : template
  /// Comparison of the second member of a std::pair<int,double>
  struct pairComparaison {
    bool operator()(const std::pair<int, double> &left,
                    const std::pair<int, double> &right) const {
      return (left.second < right.second);
    }
  };

  class MandatoryCriticalPoints : virtual public Debug {

  public:
    enum class PointType : unsigned char {
      Minimum = 0,
      JoinSaddle = 1,
      SplitSaddle = 2,
      Maximum = 3
    };

    enum class TreeType { JoinTree, SplitTree };
    MandatoryCriticalPoints();

  public:
    inline int buildJoinTreePlanarLayout() {
      return computePlanarLayout(
        TreeType::JoinTree, mdtJoinTree_, mdtJoinTreePointType_,
        mdtJoinTreePointLowInterval_, mdtJoinTreePointUpInterval_,
        mdtJoinTreePointXCoord_, mdtJoinTreePointYCoord_);
    }

    inline int buildSplitTreePlanarLayout() {
      return computePlanarLayout(
        TreeType::SplitTree, mdtSplitTree_, mdtSplitTreePointType_,
        mdtSplitTreePointLowInterval_, mdtSplitTreePointUpInterval_,
        mdtSplitTreePointXCoord_, mdtSplitTreePointYCoord_);
    }

    /// Process the construction of the 4 trees :
    /// \li Join tree of the upper bound.
    /// \li Join tree of the lower bound.
    /// \li Split tree of the upper bound.
    /// \li Split tree of the lower bound.
    /// \pre To build these trees, the following must have been called :
    /// setVertexNumber(), fillVertexScalars(), setVertexPosition(),
    /// setTriangulation() and setSoSoffsets().
    template <typename triangulationType>
    int buildSubTrees(const triangulationType &triangulation);

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    template <class dataType, typename triangulationType>
    int execute(const triangulationType &triangulation);

    template <class dataType>
    inline int fillVertexScalars(void *upperData, void *lowerData) {
      dataType *upperBoundField = (dataType *)upperData;
      dataType *lowerBoundField = (dataType *)lowerData;
      if((int)upperVertexScalars_.size() != vertexNumber_)
        upperVertexScalars_.resize(vertexNumber_);
      if((int)lowerVertexScalars_.size() != vertexNumber_)
        lowerVertexScalars_.resize(vertexNumber_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int i = 0; i < vertexNumber_; i++) {
        upperVertexScalars_[i] = (double)upperBoundField[i];
        lowerVertexScalars_[i] = (double)lowerBoundField[i];
      }
      return 0;
    }

    int findCommonAncestorNodeId(const SubLevelSetTree *tree,
                                 const int &vertexId0,
                                 const int &vertexId1) const;

    void flush();

    inline double getGlobalMaximum() const {
      return globalMaximumValue_;
    }
    inline double getGlobalMinimum() const {
      return globalMinimumValue_;
    }

    inline const Graph *getJoinTreeGraph() {
      return &(mdtJoinTree_);
    }

    inline const std::vector<double> *getJoinTreeXLayout() {
      return &(mdtJoinTreePointXCoord_);
    }

    inline const std::vector<double> *getJoinTreeYLayout() {
      return &(mdtJoinTreePointYCoord_);
    }

    inline const Graph *getSplitTreeGraph() {
      return &(mdtSplitTree_);
    }

    inline const std::vector<double> *getSplitTreeXLayout() {
      return &(mdtSplitTreePointXCoord_);
    }

    inline const std::vector<double> *getSplitTreeYLayout() {
      return &(mdtSplitTreePointYCoord_);
    }

    inline const std::vector<int> *getMdtJoinTreePointComponentId() const {
      return &(mdtJoinTreePointComponentId_);
    }
    inline const std::vector<int> *getMdtSplitTreePointComponentId() const {
      return &(mdtSplitTreePointComponentId_);
    }
    inline const std::vector<PointType> *getMdtJoinTreePointType() const {
      return &(mdtJoinTreePointType_);
    }
    inline const std::vector<PointType> *getMdtSplitTreePointType() const {
      return &(mdtSplitTreePointType_);
    }
    inline const std::vector<double> *getMdtJoinTreePointLowInterval() const {
      return &(mdtJoinTreePointLowInterval_);
    }
    inline const std::vector<double> *getMdtSplitTreePointLowInterval() const {
      return &(mdtSplitTreePointLowInterval_);
    }
    inline const std::vector<double> *getMdtJoinTreePointUpInterval() const {
      return &(mdtJoinTreePointUpInterval_);
    }
    inline const std::vector<double> *getMdtSplitTreePointUpInterval() const {
      return &(mdtSplitTreePointUpInterval_);
    }
    inline const std::vector<int> *getMdtJoinTreeEdgeSwitchable() {
      return &(mdtJoinTreeEdgeSwitchable_);
    }
    inline const std::vector<int> *getMdtSplitTreeEdgeSwitchable() {
      return &(mdtSplitTreeEdgeSwitchable_);
    }

    inline bool areSaddlesSwitchables(const TreeType treeType,
                                      const int &firstId,
                                      const int &secondId) const {
      const double firstLower
        = (treeType == TreeType::JoinTree)
            ? lowerVertexScalars_[mandatoryJoinSaddleVertex_[firstId].first]
            : lowerVertexScalars_[mandatorySplitSaddleVertex_[firstId].first];
      const double firstUpper
        = (treeType == TreeType::JoinTree)
            ? upperVertexScalars_[mandatoryJoinSaddleVertex_[firstId].second]
            : upperVertexScalars_[mandatorySplitSaddleVertex_[firstId].second];
      const double secondLower
        = (treeType == TreeType::JoinTree)
            ? lowerVertexScalars_[mandatoryJoinSaddleVertex_[secondId].first]
            : lowerVertexScalars_[mandatorySplitSaddleVertex_[secondId].first];
      const double secondUpper
        = (treeType == TreeType::JoinTree)
            ? upperVertexScalars_[mandatoryJoinSaddleVertex_[secondId].second]
            : upperVertexScalars_[mandatorySplitSaddleVertex_[secondId].second];
      return !((secondUpper < firstLower) || (firstUpper < secondLower));
    }

    // TODO Mettre les fonctions d'output dans le cpp

    template <typename triangulationType>
    inline int computeAllJoinSaddle(const triangulationType &triangulation) {
      if(mandatoryJoinSaddleVertex_.size() > 0) {
        computeJoinSaddle(0, triangulation, true);
        for(int i = 1; i < (int)mandatoryJoinSaddleVertex_.size(); i++)
          computeJoinSaddle(i, triangulation, false);
      } else {
        int *output = outputMandatoryJoinSaddle_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    int computeAllMaxima() {
      if(mandatoryMaximumVertex_.size() > 0) {
        computeMaximum(0, true, false);
        for(int i = 0; i < (int)mandatoryMaximumVertex_.size(); i++)
          computeMaximum(i, false, false);
      } else {
        int *output = outputMandatoryMaximum_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    int computeAllMinima() {
      if(mandatoryMinimumVertex_.size() > 0) {
        computeMinimum(0, true, false);
        for(int i = 0; i < (int)mandatoryMinimumVertex_.size(); i++)
          computeMinimum(i, false, false);
      } else {
        int *output = outputMandatoryMinimum_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    template <typename triangulationType>
    inline int computeAllSplitSaddle(const triangulationType &triangulation) {
      if(mandatorySplitSaddleVertex_.size() > 0) {
        computeSplitSaddle(0, triangulation, true);
        for(int i = 1; i < (int)mandatorySplitSaddleVertex_.size(); i++)
          computeSplitSaddle(i, triangulation, false);
      } else {
        int *output = outputMandatorySplitSaddle_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    template <typename triangulationType>
    inline int computeJoinSaddle(const int &id,
                                 const triangulationType &triangulation,
                                 const bool &reset = true) {
      int *output = outputMandatoryJoinSaddle_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      if(id < (int)mandatoryJoinSaddleVertex_.size()) {
        if(!isMdtJoinSaddleSimplified_[id]) {
          if(mandatoryJoinSaddleComponentVertices_[id].empty())
            computeSaddleComponent(
              id, PointType::JoinSaddle, mandatoryJoinSaddleVertex_,
              lowerVertexScalars_, upperVertexScalars_,
              mandatoryJoinSaddleComponentVertices_[id], triangulation);
          for(int i = 0;
              i < (int)mandatoryJoinSaddleComponentVertices_[id].size(); i++) {
            output[mandatoryJoinSaddleComponentVertices_[id][i]] += 1;
          }
        }
      }
      return 0;
    }

    int computeMaximum(const int &id,
                       const bool &reset = true,
                       const bool &ttkNotUsed(parallel) = true) {
      int *output = outputMandatoryMaximum_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++)
          output[i] = -1;
      if(id < (int)mandatoryMaximumVertex_.size()) {
        if(!isMdtMaximumSimplified_[id]) {
          if(mandatoryMaximumComponentVertices_[id].empty())
            computeExtremumComponent(
              PointType::Maximum, upperSplitTree_, mandatoryMaximumVertex_[id],
              lowerVertexScalars_, mandatoryMaximumComponentVertices_[id]);
          for(int i = 0; i < (int)mandatoryMaximumComponentVertices_[id].size();
              i++)
            output[mandatoryMaximumComponentVertices_[id][i]] = id;
        }
      }
      return 0;
    }

    int computeMinimum(const int &id,
                       const bool &reset = true,
                       const bool &ttkNotUsed(parallel) = true) {
      int *output = outputMandatoryMinimum_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++)
          output[i] = -1;
      if(id < (int)mandatoryMinimumVertex_.size()) {
        if(!isMdtMinimumSimplified_[id]) {
          if(mandatoryMinimumComponentVertices_[id].empty())
            computeExtremumComponent(
              PointType::Minimum, lowerJoinTree_, mandatoryMinimumVertex_[id],
              upperVertexScalars_, mandatoryMinimumComponentVertices_[id]);
          for(int i = 0; i < (int)mandatoryMinimumComponentVertices_[id].size();
              i++)
            output[mandatoryMinimumComponentVertices_[id][i]] = id;
        }
      }
      return 0;
    }

    template <typename triangulationType>
    inline int computeSplitSaddle(const int &id,
                                  const triangulationType &triangulation,
                                  const bool &reset = true) {
      int *output = outputMandatorySplitSaddle_;
      if(reset) {
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      if(id < (int)mandatorySplitSaddleVertex_.size()) {
        if(!isMdtSplitSaddleSimplified_[id]) {
          if(mandatorySplitSaddleComponentVertices_[id].empty())
            computeSaddleComponent(
              id, PointType::SplitSaddle, mandatorySplitSaddleVertex_,
              lowerVertexScalars_, upperVertexScalars_,
              mandatorySplitSaddleComponentVertices_[id], triangulation);
          for(int i = 0;
              i < (int)mandatorySplitSaddleComponentVertices_[id].size(); i++) {
            output[mandatorySplitSaddleComponentVertices_[id][i]] += 1;
          }
        }
      }
      return 0;
    }

    inline int setDebugLevel(const int &debugLevel) {
      Debug::setDebugLevel(debugLevel);
      upperJoinTree_.setDebugLevel(debugLevel);
      lowerJoinTree_.setDebugLevel(debugLevel);
      upperSplitTree_.setDebugLevel(debugLevel);
      lowerSplitTree_.setDebugLevel(debugLevel);
      return 0;
    }

    inline int setLowerBoundFieldPointer(void *data) {
      inputLowerBoundField_ = data;
      return 0;
    }

    inline int setOutputJoinSaddleDataPointer(void *data) {
      outputMandatoryJoinSaddle_ = static_cast<int *>(data);
      return 0;
    }

    /// Pass a pointer to an output array representing a scalar field.
    /// The scalars are the ids of the mandatory maximum components.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber()
    inline int setOutputMaximumDataPointer(void *data) {
      outputMandatoryMaximum_ = static_cast<int *>(data);
      return 0;
    }

    /// Pass a pointer to an output array representing a scalar field.
    /// The scalars are the ids of the mandatory minimum components.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber().
    inline int setOutputMinimumDataPointer(void *data) {
      outputMandatoryMinimum_ = static_cast<int *>(data);
      return 0;
    }

    inline int setOutputSplitSaddleDataPointer(void *data) {
      outputMandatorySplitSaddle_ = static_cast<int *>(data);
      return 0;
    }

    inline int setUpperBoundFieldPointer(void *data) {
      inputUpperBoundField_ = data;
      return 0;
    }

    inline int setSimplificationThreshold(double normalizedThreshold) {
      normalizedThreshold_ = normalizedThreshold;
      return 0;
    }

    inline int setSoSoffsets(int *offsets = NULL) {
      if((int)vertexSoSoffsets_.size() != vertexNumber_)
        vertexSoSoffsets_.resize(vertexNumber_);
      if(offsets) {
        for(int i = 0; i < vertexNumber_; i++) {
          vertexSoSoffsets_[i] = offsets[i];
        }
      } else {
        for(int i = 0; i < vertexNumber_; i++) {
          vertexSoSoffsets_[i] = i;
        }
      }
      return 0;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
      }
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    inline int setVertexNumber(const int &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    /// Set the position (x,y,z) of the i-th point
    /// \param i Index of the vertex
    /// \param point Position (x,y,z) (buffer of 3 doubles)
    /// \return Returns 0 upon success, negative values otherwise.
    inline int setVertexPosition(const int &i, const double point[3]) {
      if((int)vertexPositions_.size() != vertexNumber_)
        vertexPositions_.resize(vertexNumber_, std::vector<double>(3));
      if(i < vertexNumber_) {
        if(vertexPositions_[i].size() != 3)
          vertexPositions_[i].resize(3);
        vertexPositions_[i][0] = point[0];
        vertexPositions_[i][1] = point[1];
        vertexPositions_[i][2] = point[2];
        return 0;
      }
      return -1;
    }

    inline int simplifyJoinTree() {
      if(simplify(normalizedThreshold_, TreeType::JoinTree,
                  mdtMinJoinSaddlePair_, mergedMinimaId_,
                  mandatoryMinimumVertex_.size(), isMdtMinimumSimplified_,
                  isMdtJoinSaddleSimplified_, mdtMinimumParentSaddleId_,
                  mdtJoinSaddleParentSaddleId_)
         == 0) {
        if(buildMandatoryTree(
             TreeType::JoinTree, mdtJoinTree_, mdtJoinTreePointComponentId_,
             mdtJoinTreePointType_, mdtJoinTreePointLowInterval_,
             mdtJoinTreePointUpInterval_, mdtJoinTreeEdgeSwitchable_,
             mdtMinimumParentSaddleId_, mdtJoinSaddleParentSaddleId_,
             isMdtMinimumSimplified_, isMdtJoinSaddleSimplified_,
             mandatoryMinimumInterval_, mandatoryJoinSaddleVertex_,
             mandatoryMinimumVertex_.size(), mandatoryJoinSaddleVertex_.size(),
             PointType::Minimum, PointType::JoinSaddle, PointType::Maximum,
             getGlobalMaximum())
           == 0) {
          return 0;
        } else
          return -1;
      } else
        return -1;
    }

    inline int simplifySplitTree() {
      if(simplify(normalizedThreshold_, TreeType::SplitTree,
                  mdtMaxSplitSaddlePair_, mergedMaximaId_,
                  mandatoryMaximumVertex_.size(), isMdtMaximumSimplified_,
                  isMdtSplitSaddleSimplified_, mdtMaximumParentSaddleId_,
                  mdtSplitSaddleParentSaddleId_)
         == 0) {
        if(buildMandatoryTree(
             TreeType::SplitTree, mdtSplitTree_, mdtSplitTreePointComponentId_,
             mdtSplitTreePointType_, mdtSplitTreePointLowInterval_,
             mdtSplitTreePointUpInterval_, mdtSplitTreeEdgeSwitchable_,
             mdtMaximumParentSaddleId_, mdtSplitSaddleParentSaddleId_,
             isMdtMaximumSimplified_, isMdtSplitSaddleSimplified_,
             mandatoryMaximumInterval_, mandatorySplitSaddleVertex_,
             mandatoryMaximumVertex_.size(), mandatorySplitSaddleVertex_.size(),
             PointType::Maximum, PointType::SplitSaddle, PointType::Minimum,
             getGlobalMinimum())
           == 0) {
          return 0;
        } else
          return -1;
      } else
        return -1;
    }

  protected:
    int buildMandatoryTree(
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
      const double globalOtherExtremumValue) const;

    /// TODO : Replace SubLevelSetTrees by scalar fields for vertex value
    int
      buildPairs(const TreeType treeType,
                 const std::vector<std::pair<int, int>> &saddleList,
                 const std::vector<std::vector<int>> &mergedExtrema,
                 const std::vector<std::pair<double, double>> &extremumInterval,
                 SubLevelSetTree &lowerTree,
                 SubLevelSetTree &upperTree,
                 std::vector<std::pair<std::pair<int, int>, double>>
                   &extremaSaddlePair) const;

    int computePlanarLayout(const TreeType &treeType,
                            const Graph &mdtTree,
                            const std::vector<PointType> &mdtTreePointType,
                            const std::vector<double> &mdtTreePointLowInterval,
                            const std::vector<double> &mdtTreePointUpInterval,
                            std::vector<double> &xCoord,
                            std::vector<double> &yCoord) const;

    int computeExtremumComponent(const PointType &pointType,
                                 const SubLevelSetTree &tree,
                                 const int seedVertexId,
                                 const std::vector<double> &vertexScalars,
                                 std::vector<int> &componentVertexList) const;

    template <typename triangulationType>
    int computeSaddleComponent(
      const int componentId,
      const PointType &pointType,
      const std::vector<std::pair<int, int>> &mandatorySaddleVertex,
      const std::vector<double> &lowVertexScalars,
      const std::vector<double> &upVertexInterval,
      std::vector<int> &componentVertexList,
      const triangulationType &triangulation) const;

    int enumerateMandatoryExtrema(
      const PointType pointType,
      SubLevelSetTree &firstTree,
      SubLevelSetTree &secondTree,
      std::vector<int> &mandatoryExtremum,
      std::vector<std::pair<double, double>> &criticalInterval) const;

    /// TODO : Multiplicity
    int enumerateMandatorySaddles(
      const PointType pointType,
      SubLevelSetTree &lowerTree,
      SubLevelSetTree &upperTree,
      const std::vector<int> &mandatoryExtremumVertex,
      std::vector<std::pair<int, int>> &mandatorySaddleVertex,
      std::vector<std::vector<int>> &mandatoryMergedExtrema);

    int getSubTreeRootSuperArcId(const SubLevelSetTree *tree,
                                 const int &startingSuperArcId,
                                 const double &targetValue) const;

    void getSubTreeSuperArcIds(const SubLevelSetTree *tree,
                               const int &rootSuperArcId,
                               std::vector<int> &subTreeSuperArcId) const;

    inline int getVertexSuperArcId(const int &vertexId,
                                   const SubLevelSetTree *tree) const {

      int superArcId = tree->getVertexSuperArcId(vertexId);
      // If superArcId == -1, it may be a leaf so look for the nearest super arc
      if(superArcId == -1) {
        int nodeId = tree->getVertexNodeId(vertexId);
        if(tree->getNode(nodeId)->getNumberOfUpSuperArcs()) {
          superArcId = tree->getNode(nodeId)->getUpSuperArcId(0);
        } else if(tree->getNode(nodeId)->getNumberOfDownSuperArcs()) {
          superArcId = tree->getNode(nodeId)->getDownSuperArcId(0);
        }
      }
      return superArcId;
    }

    int simplify(const double normalizedThreshold,
                 const TreeType treeType,
                 const std::vector<std::pair<std::pair<int, int>, double>>
                   &extremaSaddlePair,
                 const std::vector<std::vector<int>> &mergedExtrema,
                 const int numberOfExtrema,
                 std::vector<bool> &extremumSimplified,
                 std::vector<bool> &saddleSimplified,
                 std::vector<int> &extremumParentSaddle,
                 std::vector<int> &saddleParentSaddle) const;

  protected:
    /// Void pointer to the input upper bound scalar field.
    void *inputUpperBoundField_{};
    /// Void pointer to the input lower bound scalar field.
    void *inputLowerBoundField_{};
    /// Void pointer to the output mandatory minima components.
    int *outputMandatoryMinimum_{};
    /// Void pointer to the output mandatory join saddles components.
    int *outputMandatoryJoinSaddle_{};
    /// Void pointer to the output mandatory split saddles components.
    int *outputMandatorySplitSaddle_{};
    /// Void pointer to the output mandatory maxima components.
    int *outputMandatoryMaximum_{};
    /// Number of vertices
    int vertexNumber_{};
    /// Position (x,y,z) of each vertex
    std::vector<std::vector<double>> vertexPositions_{};
    /// Offsets
    std::vector<int> vertexSoSoffsets_{};
    /// Copy of the input upper scalar field converted in double.
    std::vector<double> upperVertexScalars_{};
    /// Copy of the input lower scalar field converted in double.
    std::vector<double> lowerVertexScalars_{};
    /// Join tree of the upper bound scalar field.
    SubLevelSetTree upperJoinTree_{};
    /// Join tree of the lower bound scalar field.
    SubLevelSetTree lowerJoinTree_{};
    /// Split tree of the upper bound scalar field.
    SubLevelSetTree upperSplitTree_{};
    /// Split tree of the lower bound scalar field.
    SubLevelSetTree lowerSplitTree_{};
    /// List of vertex id of the minima in the upper bound scalar field.
    std::vector<int> upperMinimumList_{};
    /// List of vertex id of the minima in the lower bound scalar field.
    std::vector<int> lowerMinimumList_{};
    /// List of vertex id of the maxima in the upper bound scalar field.
    std::vector<int> upperMaximumList_{};
    /// List of vertex id of the maxima in the lower bound scalar field.
    std::vector<int> lowerMaximumList_{};
    /// Mandatory vertex for each minimum component.
    std::vector<int> mandatoryMinimumVertex_{};
    /// Mandatory vertex for each maximum component.
    std::vector<int> mandatoryMaximumVertex_{};
    /// Critical interval for each minimum component
    std::vector<std::pair<double, double>> mandatoryMinimumInterval_{};
    /// Critical interval for each maximum component
    std::vector<std::pair<double, double>> mandatoryMaximumInterval_{};
    /// Pair of mandatory vertices for each join saddle component.
    std::vector<std::pair<int, int>> mandatoryJoinSaddleVertex_{};
    /// Pair of mandatory vertices for each split saddle component.
    std::vector<std::pair<int, int>> mandatorySplitSaddleVertex_{};
    /// List of ids of the mandatory minima merged for each mandatory join
    /// saddle.
    std::vector<std::vector<int>> mergedMaximaId_{};
    /// List of ids of the mandatory maxima merged for each mandatory split
    /// saddle.
    std::vector<std::vector<int>> mergedMinimaId_{};
    /// Pairs ( (M,S), d(M,S) ) Of minima and join saddles
    std::vector<std::pair<std::pair<int, int>, double>> mdtMinJoinSaddlePair_{};
    /// Pairs ( (M,S), d(M,S) ) Of maxima and split saddles
    std::vector<std::pair<std::pair<int, int>, double>>
      mdtMaxSplitSaddlePair_{};
    /// Value of the simplification threshold.
    double normalizedThreshold_{0.0};
    /// Flags indicating if the mandatory minimum component have been
    /// simplified.
    std::vector<bool> isMdtMinimumSimplified_;
    /// Flags indicating if the mandatory join saddle component have been
    /// simplified.
    std::vector<bool> isMdtJoinSaddleSimplified_{};
    /// Flags indicating if the mandatory split saddle component have been
    /// simplified.
    std::vector<bool> isMdtSplitSaddleSimplified_{};
    /// Flags indicating if the mandatory maximum component have been
    /// simplified.
    std::vector<bool> isMdtMaximumSimplified_{};

    // Graph
    std::vector<int> mdtMinimumParentSaddleId_{};
    std::vector<int> mdtJoinSaddleParentSaddleId_{};
    std::vector<int> mdtSplitSaddleParentSaddleId_{};
    std::vector<int> mdtMaximumParentSaddleId_{};

    Graph mdtJoinTree_{};
    Graph mdtSplitTree_{};

    std::vector<int> mdtJoinTreePointComponentId_{};
    std::vector<int> mdtSplitTreePointComponentId_{};
    std::vector<PointType> mdtJoinTreePointType_{};
    std::vector<PointType> mdtSplitTreePointType_{};
    std::vector<double> mdtJoinTreePointLowInterval_{};
    std::vector<double> mdtSplitTreePointLowInterval_{};
    std::vector<double> mdtJoinTreePointUpInterval_{};
    std::vector<double> mdtSplitTreePointUpInterval_{};

    std::vector<int> mdtJoinTreeEdgeSwitchable_{};
    std::vector<int> mdtSplitTreeEdgeSwitchable_{};

    std::vector<double> mdtJoinTreePointXCoord_{};
    std::vector<double> mdtSplitTreePointXCoord_{};

    std::vector<double> mdtJoinTreePointYCoord_{};
    std::vector<double> mdtSplitTreePointYCoord_{};

    double globalMinimumValue_{};
    double globalMaximumValue_{};

    /// List of the vertices forming each of the mandatory maximum components.
    std::vector<std::vector<int>> mandatoryMaximumComponentVertices_{};
    /// List of the vertices forming each of the mandatory split saddle
    /// components.
    std::vector<std::vector<int>> mandatoryMinimumComponentVertices_{};
    /// List of the vertices forming each of the mandatory join saddle
    /// components.
    std::vector<std::vector<int>> mandatoryJoinSaddleComponentVertices_{};
    /// List of the vertices forming each of the mandatory split saddle
    /// components.
    std::vector<std::vector<int>> mandatorySplitSaddleComponentVertices_{};
  };
} // namespace ttk

// if the package is a pure template class, uncomment the following line
// #include                  <MandatoryCriticalPoints.cpp>

// template functions
template <class dataType, typename triangulationType>
int ttk::MandatoryCriticalPoints::execute(
  const triangulationType &triangulation) {

  Timer t;

// Check the consistency of the variables
// TODO Move the output checks in the right functions
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputUpperBoundField_)
    return -1;
  if(!inputLowerBoundField_)
    return -2;
  if(!outputMandatoryMinimum_)
    return -3;
  if(!outputMandatoryJoinSaddle_)
    return -4;
  if(!outputMandatorySplitSaddle_)
    return -5;
  if(!outputMandatoryMaximum_)
    return -6;
  if(!vertexNumber_)
    return -7;
  if(!vertexPositions_.size())
    return -8;
#endif

  // Init the input
  fillVertexScalars<dataType>(inputUpperBoundField_, inputLowerBoundField_);

  // Build the join trees and split trees
  buildSubTrees(triangulation);

// Compute mandatory extrema
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      enumerateMandatoryExtrema(PointType::Minimum, upperJoinTree_,
                                lowerJoinTree_, mandatoryMinimumVertex_,
                                mandatoryMinimumInterval_);
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      enumerateMandatoryExtrema(PointType::Maximum, lowerSplitTree_,
                                upperSplitTree_, mandatoryMaximumVertex_,
                                mandatoryMaximumInterval_);
    }
  }

  // Compute mandatory saddles
  enumerateMandatorySaddles(PointType::JoinSaddle, lowerJoinTree_,
                            upperJoinTree_, mandatoryMinimumVertex_,
                            mandatoryJoinSaddleVertex_, mergedMinimaId_);
  enumerateMandatorySaddles(PointType::SplitSaddle, lowerSplitTree_,
                            upperSplitTree_, mandatoryMaximumVertex_,
                            mandatorySplitSaddleVertex_, mergedMaximaId_);

  // Build pairs of <extremum,saddle>
  buildPairs(TreeType::JoinTree, mandatoryJoinSaddleVertex_, mergedMinimaId_,
             mandatoryMinimumInterval_, lowerJoinTree_, upperJoinTree_,
             mdtMinJoinSaddlePair_);
  buildPairs(TreeType::SplitTree, mandatorySplitSaddleVertex_, mergedMaximaId_,
             mandatoryMaximumInterval_, lowerSplitTree_, upperSplitTree_,
             mdtMaxSplitSaddlePair_);

  // Simplify pairs
  simplify(normalizedThreshold_, TreeType::JoinTree, mdtMinJoinSaddlePair_,
           mergedMinimaId_, mandatoryMinimumVertex_.size(),
           isMdtMinimumSimplified_, isMdtJoinSaddleSimplified_,
           mdtMinimumParentSaddleId_, mdtJoinSaddleParentSaddleId_);
  simplify(normalizedThreshold_, TreeType::SplitTree, mdtMaxSplitSaddlePair_,
           mergedMaximaId_, mandatoryMaximumVertex_.size(),
           isMdtMaximumSimplified_, isMdtSplitSaddleSimplified_,
           mdtMaximumParentSaddleId_, mdtSplitSaddleParentSaddleId_);

  // Build the mandatory trees
  buildMandatoryTree(TreeType::JoinTree, mdtJoinTree_,
                     mdtJoinTreePointComponentId_, mdtJoinTreePointType_,
                     mdtJoinTreePointLowInterval_, mdtJoinTreePointUpInterval_,
                     mdtJoinTreeEdgeSwitchable_, mdtMinimumParentSaddleId_,
                     mdtJoinSaddleParentSaddleId_, isMdtMinimumSimplified_,
                     isMdtJoinSaddleSimplified_, mandatoryMinimumInterval_,
                     mandatoryJoinSaddleVertex_, mandatoryMinimumVertex_.size(),
                     mandatoryJoinSaddleVertex_.size(), PointType::Minimum,
                     PointType::JoinSaddle, PointType::Maximum,
                     getGlobalMaximum());
  buildMandatoryTree(
    TreeType::SplitTree, mdtSplitTree_, mdtSplitTreePointComponentId_,
    mdtSplitTreePointType_, mdtSplitTreePointLowInterval_,
    mdtSplitTreePointUpInterval_, mdtSplitTreeEdgeSwitchable_,
    mdtMaximumParentSaddleId_, mdtSplitSaddleParentSaddleId_,
    isMdtMaximumSimplified_, isMdtSplitSaddleSimplified_,
    mandatoryMaximumInterval_, mandatorySplitSaddleVertex_,
    mandatoryMaximumVertex_.size(), mandatorySplitSaddleVertex_.size(),
    PointType::Maximum, PointType::SplitSaddle, PointType::Minimum,
    getGlobalMinimum());

  // Compute the planar layout for the output trees
  computePlanarLayout(TreeType::JoinTree, mdtJoinTree_, mdtJoinTreePointType_,
                      mdtJoinTreePointLowInterval_, mdtJoinTreePointUpInterval_,
                      mdtJoinTreePointXCoord_, mdtJoinTreePointYCoord_);
  computePlanarLayout(TreeType::SplitTree, mdtSplitTree_,
                      mdtSplitTreePointType_, mdtSplitTreePointLowInterval_,
                      mdtSplitTreePointUpInterval_, mdtSplitTreePointXCoord_,
                      mdtSplitTreePointYCoord_);

  // Clear outputs
  mandatoryMinimumComponentVertices_.resize(mandatoryMinimumVertex_.size());
  fill(mandatoryMinimumComponentVertices_.begin(),
       mandatoryMinimumComponentVertices_.end(), std::vector<int>());
  mandatoryJoinSaddleComponentVertices_.resize(
    mandatoryJoinSaddleVertex_.size());
  fill(mandatoryJoinSaddleComponentVertices_.begin(),
       mandatoryJoinSaddleComponentVertices_.end(), std::vector<int>());
  mandatorySplitSaddleComponentVertices_.resize(
    mandatorySplitSaddleVertex_.size());
  fill(mandatorySplitSaddleComponentVertices_.begin(),
       mandatorySplitSaddleComponentVertices_.end(), std::vector<int>());
  mandatoryMaximumComponentVertices_.resize(mandatoryMaximumVertex_.size());
  fill(mandatoryMaximumComponentVertices_.begin(),
       mandatoryMaximumComponentVertices_.end(), std::vector<int>());

  this->printMsg(
    "Data-set (" + std::to_string(vertexNumber_) + " points) processed", 1.0,
    t.getElapsedTime(), this->threadNumber_);
  return 0;
}

template <typename triangulationType>
int ttk::MandatoryCriticalPoints::buildSubTrees(
  const triangulationType &triangulation) {

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
  if(triangulation.getNumberOfVertices() != vertexNumber_)
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
    SimplexId neighborNumber = triangulation.getVertexNeighborNumber(i);
    for(SimplexId j = 0; j < neighborNumber; j++) {
      SimplexId neighborId;
      triangulation.getVertexNeighbor(i, j, neighborId);
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
    {
      upperJoinTree_.setNumberOfVertices(vertexNumber_);
      upperJoinTree_.setVertexScalars(&upperVertexScalars_);
      upperJoinTree_.setVertexPositions(&vertexPositions_);
      upperJoinTree_.setTriangulation(&triangulation);
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
      lowerJoinTree_.setTriangulation(&triangulation);
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
      upperSplitTree_.setTriangulation(&triangulation);
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
      lowerSplitTree_.setTriangulation(&triangulation);
      lowerSplitTree_.setVertexSoSoffsets(&vertexSoSoffsets_);
      lowerSplitTree_.buildExtremumList(lowerMaximumList_, false);
      lowerSplitTree_.build();
    }
  }
  this->printMsg("4 SubLevelSetTrees computed", 1.0, t.getElapsedTime(),
                 this->threadNumber_);
  return 0;
}

template <typename triangulationType>
int ttk::MandatoryCriticalPoints::computeSaddleComponent(
  const int componentId,
  const PointType &pointType,
  const std::vector<std::pair<int, int>> &mandatorySaddleVertex,
  const std::vector<double> &lowVertexScalars,
  const std::vector<double> &upVertexScalars,
  std::vector<int> &componentVertexList,
  const triangulationType &triangulation) const {

  const int seedVertexId = mandatorySaddleVertex[componentId].first;
  const double lowInterval = lowVertexScalars[seedVertexId];
  const double upInterval
    = upVertexScalars[mandatorySaddleVertex[componentId].second];

  componentVertexList.clear();

  std::vector<bool> isVisited(vertexNumber_, false);
  std::queue<SimplexId> idQueue;
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
        componentVertexList.push_back(vertexId);
        // Neighbors
        SimplexId neighborNumber
          = triangulation.getVertexNeighborNumber(vertexId);
        for(SimplexId i = 0; i < neighborNumber; i++) {
          SimplexId neighborVertexId;
          triangulation.getVertexNeighbor(vertexId, i, neighborVertexId);
          idQueue.push(neighborVertexId);
        }
      }
    }
  }
  return 0;
}
