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

#ifndef _MANDATORYCRITICALPOINTS_H
#define _MANDATORYCRITICALPOINTS_H

// base code includes
#include <ContourTree.h>
#include <LowestCommonAncestor.h>
#include <Triangulation.h>
#include <Wrapper.h>

// std includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <queue>
#include <utility>
#include <vector>

namespace ttk {

  class Graph : public Debug {
  public:
    class Vertex {
      friend class Graph;

    public:
      Vertex() {
      }
      ~Vertex() {
      }
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
      std::vector<int> edgeIdx_;
    };
    class Edge {
      friend class Graph;

    public:
      Edge(const int &start, const int &end) {
        vertexIdx_.first = start;
        vertexIdx_.second = end;
      }
      ~Edge() {
      }
      const std::pair<int, int> &getVertexIdx() const {
        return vertexIdx_;
      }

    protected:
      std::pair<int, int> vertexIdx_;
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
    const Vertex *getVertex(const int &idx) const {
      if(idx < (int)vertexList_.size())
        return &(vertexList_[idx]);
      else
        return NULL;
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
    const Edge *getEdge(const int &idx) const {
      if(idx < (int)edgeList_.size())
        return &(edgeList_[idx]);
      else
        return NULL;
    }
    inline void clear() {
      vertexList_.clear();
      edgeList_.clear();
    }

  protected:
    std::vector<Vertex> vertexList_;
    std::vector<Edge> edgeList_;
  };

  /// Comparison between critical point pairs ( (Extremum,Saddle), dist(M,S) )
  struct criticalPointPairComparaison {
    bool operator()(const std::pair<std::pair<int, int>, double> &left,
                    const std::pair<std::pair<int, int>, double> &right) {
      return (left.second < right.second);
    }
  };

  /// Comparison between mandatory saddles (Saddle id, Number of merged extrema)
  struct mandatorySaddleComparaison {
    bool operator()(const std::pair<int, int> &left,
                    const std::pair<int, int> &right) {
      return (left.second < right.second);
    }
  };

  // TODO : template
  /// Comparison of the second member of a std::pair<int,double>
  struct pairComparaison {
    bool operator()(const std::pair<int, double> &left,
                    const std::pair<int, double> &right) {
      return (left.second < right.second);
    }
  };

  class MandatoryCriticalPoints : public Debug {

  public:
    enum PointType : unsigned char {
      Minimum = 0,
      JoinSaddle = 1,
      SplitSaddle = 2,
      Maximum = 3
    };

    enum TreeType { JoinTree, SplitTree };

    MandatoryCriticalPoints();
    ~MandatoryCriticalPoints();

  public:
    inline int buildJoinTreePlanarLayout() {
      return computePlanarLayout(TreeType::JoinTree);
    }

    inline int buildSplitTreePlanarLayout() {
      return computePlanarLayout(TreeType::SplitTree);
    }

    /// Process the construction of the 4 trees :
    /// \li Join tree of the upper bound.
    /// \li Join tree of the lower bound.
    /// \li Split tree of the upper bound.
    /// \li Split tree of the lower bound.
    /// \pre To build these trees, the following must have been called :
    /// setVertexNumber(), fillVertexScalars(), setVertexPosition(),
    /// setTriangulation() and setSoSoffsets().
    int buildSubTrees();

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    template <class dataType>
    int execute();

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

    inline int outputAllJoinSaddle() {
      if(mandatoryJoinSaddleVertex_.size() > 0) {
        outputJoinSaddle(0, true);
        for(int i = 1; i < (int)mandatoryJoinSaddleVertex_.size(); i++)
          outputJoinSaddle(i, false);
      } else {
        int *output = (int *)outputMandatoryJoinSaddle_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    int outputAllMaxima() {
      if(mandatoryMaximumVertex_.size() > 0) {
        outputMaximum(0, true, false);
        for(int i = 0; i < (int)mandatoryMaximumVertex_.size(); i++)
          outputMaximum(i, false, false);
      } else {
        int *output = (int *)outputMandatoryMaximum_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    int outputAllMinima() {
      if(mandatoryMinimumVertex_.size() > 0) {
        outputMinimum(0, true, false);
        for(int i = 0; i < (int)mandatoryMinimumVertex_.size(); i++)
          outputMinimum(i, false, false);
      } else {
        int *output = (int *)outputMandatoryMinimum_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    inline int outputAllSplitSaddle() {
      if(mandatorySplitSaddleVertex_.size() > 0) {
        outputSplitSaddle(0, true);
        for(int i = 1; i < (int)mandatorySplitSaddleVertex_.size(); i++)
          outputSplitSaddle(i, false);
      } else {
        int *output = (int *)outputMandatorySplitSaddle_;
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      return 0;
    }

    inline int outputJoinSaddle(const int &id, const bool &reset = true) {
      int *output = (int *)outputMandatoryJoinSaddle_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      if(id < (int)mandatoryJoinSaddleVertex_.size()) {
        if(!isMdtJoinSaddleSimplified_[id]) {
          if(mandatoryJoinSaddleComponentVertices_[id].empty())
            computeSaddleComponent(id, PointType::JoinSaddle);
          for(int i = 0;
              i < (int)mandatoryJoinSaddleComponentVertices_[id].size(); i++) {
            output[mandatoryJoinSaddleComponentVertices_[id][i]] += 1;
          }
        }
      }
      return 0;
    }

    int outputMaximum(const int &id,
                      const bool &reset = true,
                      const bool &parallel = true) {
      int *output = (int *)outputMandatoryMaximum_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++)
          output[i] = -1;
      if(id < (int)mandatoryMaximumVertex_.size()) {
        if(!isMdtMaximumSimplified_[id]) {
          if(mandatoryMaximumComponentVertices_[id].empty())
            computeExtremumComponent(id, PointType::Maximum);
          for(int i = 0; i < (int)mandatoryMaximumComponentVertices_[id].size();
              i++)
            output[mandatoryMaximumComponentVertices_[id][i]] = id;
        }
      }
      return 0;
    }

    int outputMinimum(const int &id,
                      const bool &reset = true,
                      const bool &parallel = true) {
      int *output = (int *)outputMandatoryMinimum_;
      if(reset)
        for(int i = 0; i < vertexNumber_; i++)
          output[i] = -1;
      if(id < (int)mandatoryMinimumVertex_.size()) {
        if(!isMdtMinimumSimplified_[id]) {
          if(mandatoryMinimumComponentVertices_[id].empty())
            computeExtremumComponent(id, PointType::Minimum);
          for(int i = 0; i < (int)mandatoryMinimumComponentVertices_[id].size();
              i++)
            output[mandatoryMinimumComponentVertices_[id][i]] = id;
        }
      }
      return 0;
    }

    inline int outputSplitSaddle(const int &id, const bool &reset = true) {
      int *output = (int *)outputMandatorySplitSaddle_;
      if(reset) {
        for(int i = 0; i < vertexNumber_; i++) {
          output[i] = -1;
        }
      }
      if(id < (int)mandatorySplitSaddleVertex_.size()) {
        if(!isMdtSplitSaddleSimplified_[id]) {
          if(mandatorySplitSaddleComponentVertices_[id].empty())
            computeSaddleComponent(id, PointType::SplitSaddle);
          for(int i = 0;
              i < (int)mandatorySplitSaddleComponentVertices_[id].size(); i++) {
            output[mandatorySplitSaddleComponentVertices_[id][i]] += 1;
          }
        }
      }
      return 0;
    }

    inline int setDebugLevel(const int &debugLevel) {
      upperJoinTree_.setDebugLevel(debugLevel);
      lowerJoinTree_.setDebugLevel(debugLevel);
      upperSplitTree_.setDebugLevel(debugLevel);
      lowerSplitTree_.setDebugLevel(debugLevel);
      debugLevel_ = debugLevel;
      return 0;
    }

    inline int setLowerBoundFieldPointer(void *data) {
      inputLowerBoundField_ = data;
      return 0;
    }

    inline int setOutputJoinSaddleDataPointer(void *data) {
      outputMandatoryJoinSaddle_ = data;
      return 0;
    }

    /// Pass a pointer to an output array representing a scalar field.
    /// The scalars are the ids of the mandatory maximum components.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber()
    inline int setOutputMaximumDataPointer(void *data) {
      outputMandatoryMaximum_ = data;
      return 0;
    }

    /// Pass a pointer to an output array representing a scalar field.
    /// The scalars are the ids of the mandatory minimum components.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber().
    inline int setOutputMinimumDataPointer(void *data) {
      outputMandatoryMinimum_ = data;
      return 0;
    }

    inline int setOutputSplitSaddleDataPointer(void *data) {
      outputMandatorySplitSaddle_ = data;
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

    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_) {
        triangulation_->preprocessVertexNeighbors();
      }
      return 0;
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
    /// \param point[3] Position (x,y,z)
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
      if(simplify(normalizedThreshold_, TreeType::JoinTree) == 0) {
        if(buildMandatoryTree(TreeType::JoinTree) == 0) {
          return 0;
        } else
          return -1;
      } else
        return -1;
    }

    inline int simplifySplitTree() {
      if(simplify(normalizedThreshold_, TreeType::SplitTree) == 0) {
        if(buildMandatoryTree(TreeType::SplitTree) == 0) {
          return 0;
        } else
          return -1;
      } else
        return -1;
    }

  protected:
    int buildMandatoryTree(const TreeType treeType);

    /// TODO : Replace SubLevelSetTrees by scalar fields for vertex value
    int buildPairs(const TreeType treeType);

    int computePlanarLayout(const TreeType &treeType);

    int computeExtremumComponent(const int &componentId,
                                 const PointType &pointType);

    int computeSaddleComponent(const int &componentId,
                               const PointType &pointType);

    int enumerateMandatoryExtrema(const PointType pointType);

    /// TODO : Multiplicity
    int enumerateMandatorySaddles(const PointType pointType);

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

    int simplify(const double &normalizedThreshold, const TreeType treeType);

  protected:
    /// Void pointer to the input upper bound scalar field.
    void *inputUpperBoundField_;
    /// Void pointer to the input lower bound scalar field.
    void *inputLowerBoundField_;
    /// Void pointer to the output mandatory minima components.
    void *outputMandatoryMinimum_;
    /// Void pointer to the output mandatory join saddles components.
    void *outputMandatoryJoinSaddle_;
    /// Void pointer to the output mandatory split saddles components.
    void *outputMandatorySplitSaddle_;
    /// Void pointer to the output mandatory maxima components.
    void *outputMandatoryMaximum_;
    /// Number of vertices
    int vertexNumber_;
    /// Position (x,y,z) of each vertex
    std::vector<std::vector<double>> vertexPositions_;
    /// Offsets
    std::vector<int> vertexSoSoffsets_;
    /// Copy of the input upper scalar field converted in double.
    std::vector<double> upperVertexScalars_;
    /// Copy of the input lower scalar field converted in double.
    std::vector<double> lowerVertexScalars_;
    /// Triangulation object of the input.
    Triangulation *triangulation_;
    /// Join tree of the upper bound scalar field.
    SubLevelSetTree upperJoinTree_;
    /// Join tree of the lower bound scalar field.
    SubLevelSetTree lowerJoinTree_;
    /// Split tree of the upper bound scalar field.
    SubLevelSetTree upperSplitTree_;
    /// Split tree of the lower bound scalar field.
    SubLevelSetTree lowerSplitTree_;
    /// List of vertex id of the minima in the upper bound scalar field.
    std::vector<int> upperMinimumList_;
    /// List of vertex id of the minima in the lower bound scalar field.
    std::vector<int> lowerMinimumList_;
    /// List of vertex id of the maxima in the upper bound scalar field.
    std::vector<int> upperMaximumList_;
    /// List of vertex id of the maxima in the lower bound scalar field.
    std::vector<int> lowerMaximumList_;
    /// Mandatory vertex for each minimum component.
    std::vector<int> mandatoryMinimumVertex_;
    /// Mandatory vertex for each maximum component.
    std::vector<int> mandatoryMaximumVertex_;
    /// Critical interval for each minimum component
    std::vector<std::pair<double, double>> mandatoryMinimumInterval_;
    /// Critical interval for each maximum component
    std::vector<std::pair<double, double>> mandatoryMaximumInterval_;
    /// Pair of mandatory vertices for each join saddle component.
    std::vector<std::pair<int, int>> mandatoryJoinSaddleVertex_;
    /// Pair of mandatory vertices for each split saddle component.
    std::vector<std::pair<int, int>> mandatorySplitSaddleVertex_;
    /// List of ids of the mandatory minima merged for each mandatory join
    /// saddle.
    std::vector<std::vector<int>> mergedMaximaId_;
    /// List of ids of the mandatory maxima merged for each mandatory split
    /// saddle.
    std::vector<std::vector<int>> mergedMinimaId_;
    /// Pairs ( (M,S), d(M,S) ) Of minima and join saddles
    std::vector<std::pair<std::pair<int, int>, double>> mdtMinJoinSaddlePair_;
    /// Pairs ( (M,S), d(M,S) ) Of maxima and split saddles
    std::vector<std::pair<std::pair<int, int>, double>> mdtMaxSplitSaddlePair;
    /// Value of the simplification threshold.
    double normalizedThreshold_;
    /// Flags indicating if the mandatory minimum component have been
    /// simplified.
    std::vector<bool> isMdtMinimumSimplified_;
    /// Flags indicating if the mandatory join saddle component have been
    /// simplified.
    std::vector<bool> isMdtJoinSaddleSimplified_;
    /// Flags indicating if the mandatory split saddle component have been
    /// simplified.
    std::vector<bool> isMdtSplitSaddleSimplified_;
    /// Flags indicating if the mandatory maximum component have been
    /// simplified.
    std::vector<bool> isMdtMaximumSimplified_;

    // Graph
    std::vector<int> mdtMinimumParentSaddleId_;
    std::vector<int> mdtJoinSaddleParentSaddleId_;
    std::vector<int> mdtSplitSaddleParentSaddleId_;
    std::vector<int> mdtMaximumParentSaddleId_;

    Graph mdtJoinTree_;
    Graph mdtSplitTree_;

    std::vector<int> mdtJoinTreePointComponentId_;
    std::vector<int> mdtSplitTreePointComponentId_;
    std::vector<PointType> mdtJoinTreePointType_;
    std::vector<PointType> mdtSplitTreePointType_;
    std::vector<double> mdtJoinTreePointLowInterval_;
    std::vector<double> mdtSplitTreePointLowInterval_;
    std::vector<double> mdtJoinTreePointUpInterval_;
    std::vector<double> mdtSplitTreePointUpInterval_;

    std::vector<int> mdtJoinTreeEdgeSwitchable_;
    std::vector<int> mdtSplitTreeEdgeSwitchable_;

    std::vector<double> mdtJoinTreePointXCoord_;
    std::vector<double> mdtSplitTreePointXCoord_;

    std::vector<double> mdtJoinTreePointYCoord_;
    std::vector<double> mdtSplitTreePointYCoord_;

    double globalMinimumValue_;
    double globalMaximumValue_;

    /// List of the vertices forming each of the mandatory maximum components.
    std::vector<std::vector<int>> mandatoryMaximumComponentVertices_;
    /// List of the vertices forming each of the mandatory split saddle
    /// components.
    std::vector<std::vector<int>> mandatoryMinimumComponentVertices_;
    /// List of the vertices forming each of the mandatory join saddle
    /// components.
    std::vector<std::vector<int>> mandatoryJoinSaddleComponentVertices_;
    /// List of the vertices forming each of the mandatory split saddle
    /// components.
    std::vector<std::vector<int>> mandatorySplitSaddleComponentVertices_;
  };
} // namespace ttk

// if the package is a pure template class, uncomment the following line
// #include                  <MandatoryCriticalPoints.cpp>

// template functions
template <class dataType>
int ttk::MandatoryCriticalPoints::execute() {

  Timer t;

// Check the consistency of the variables
// TODO Déplacer les vérifications des outputs dans les bonnes fonctions
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
  if(!triangulation_)
    return -9;
#endif

  // Init the input
  fillVertexScalars<dataType>(inputUpperBoundField_, inputLowerBoundField_);

  // Build the join trees and split trees
  buildSubTrees();

// Compute mandatory extrema
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { enumerateMandatoryExtrema(PointType::Minimum); }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { enumerateMandatoryExtrema(PointType::Maximum); }
  }

  // Compute mandatory saddles
  enumerateMandatorySaddles(PointType::JoinSaddle);
  enumerateMandatorySaddles(PointType::SplitSaddle);

  // Build pairs of <extremum,saddle>
  buildPairs(TreeType::JoinTree);
  buildPairs(TreeType::SplitTree);

  // Simplify pairs
  simplify(normalizedThreshold_, TreeType::JoinTree);
  simplify(normalizedThreshold_, TreeType::SplitTree);

  // Build the mandatory trees
  buildMandatoryTree(TreeType::JoinTree);
  buildMandatoryTree(TreeType::SplitTree);

  // Compute the planar layout for the output trees
  computePlanarLayout(TreeType::JoinTree);
  computePlanarLayout(TreeType::SplitTree);

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

  // Debug messages
  if(debugLevel_ > timeMsg) {
    std::stringstream msg;
    msg << "[MandatoryCriticalPoints] Data-set (" << vertexNumber_
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  return 0;
}

#endif // MANDATORYCRITICALPOINTS_H
