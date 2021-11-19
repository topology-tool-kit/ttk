/// \ingroup base
/// \class ttk::ContourTree
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date July 2011
///
/// \brief TTK processing package that computes the contour tree of scalar data
/// and more (data segmentation, topological simplification, persistence
/// diagrams, persistence curves, etc.).
///
/// \warning SimplexId (large large datasets). This class builds and runs
/// with the new triangulation API (SimplexId) but may need adjustments when
/// addressing more than integers (large datasets).
///
/// \b Related \b publication \n
/// "Computing contour trees in all dimensions" \n
/// Hamish Carr, Jack Snoeyink, Ulrike Axen \n
/// Proc. of ACM SODA 2000.
///

#pragma once

#include <Triangulation.h>
#include <UnionFind.h>

#include <vector>

namespace ttk {

  class SubLevelSetTree;
  class Node;

  class Arc : virtual public Debug {
  public:
    inline int getDownNodeId() const {
      return downNodeId_;
    }

    inline int getUpNodeId() const {
      return upNodeId_;
    }

    inline void setDownNodeId(const int &downNodeId) {
      downNodeId_ = downNodeId;
    }

    inline void setUpNodeId(const int &upNodeId) {
      upNodeId_ = upNodeId;
    }

  protected:
    int downNodeId_{-1}, upNodeId_{-1};
  };

  class SuperArc : public Arc {
  public:
    friend class SubLevelSetTree;

    inline void appendRegularNode(const int &nodeId) {
      regularNodeList_.push_back(nodeId);
    }

    inline void appendBarycenter(const std::vector<double> &barycenter) {
      barycenterList_.push_back(barycenter);
    }

    inline void appendSample(const std::vector<int> &sample) {
      sampleList_.push_back(sample);
    }

    inline int getNumberOfRegularNodes() const {
      return (int)regularNodeList_.size();
    }

    inline int getNumberOfBarycenters() const {
      return (int)barycenterList_.size();
    }

    inline int getNumberOfSamples() const {
      return (int)sampleList_.size();
    }

    inline int getRegularNodeId(const int &arcNodeId) const {
      if((arcNodeId < 0)
         || (((unsigned int)arcNodeId) >= regularNodeList_.size()))
        return -1;
      return regularNodeList_[arcNodeId];
    }

    inline void getBarycenter(const int &id,
                              std::vector<double> &barycenter) const {
      for(unsigned int k = 0; k < 3; ++k)
        barycenter[k] = barycenterList_[id][k];
    }

    inline void getBarycenter(const int &id, double barycenter[3]) const {
      for(unsigned int k = 0; k < 3; ++k)
        barycenter[k] = barycenterList_[id][k];
    }

    inline void getSample(const int &id, std::vector<int> &sample) const {
      sample = sampleList_[id];
    }

    inline void clearBarycenters() {
      barycenterList_.clear();
    }

    inline void clearSamples() {
      sampleList_.clear();
    }

    inline void clearRegularNodes() {
      regularNodeList_.clear();
    }

    inline bool isPruned() const {
      return pruned_;
    }

    void smooth(const std::vector<Node> &nodeList,
                const std::vector<std::vector<double>> *vertexPositions,
                bool order = true);

    void sortRegularNodes(const std::vector<double> *vertexScalars,
                          const std::vector<int> *vertexOffsets,
                          const std::vector<Node> *nodeList,
                          bool order = true);

  protected:
    bool pruned_{false};
    std::vector<int> regularNodeList_{};
    std::vector<std::vector<double>> barycenterList_{};
    std::vector<std::vector<int>> sampleList_{};
  };

  class Node : virtual public Debug {
  public:
    inline void addDownArcId(const int &downArcId) {
      downArcList_.push_back(downArcId);
    }

    inline void addDownSuperArcId(const int &downSuperArcId) {
      downSuperArcList_.push_back(downSuperArcId);
    }

    inline void addUpArcId(const int &upArcId) {
      upArcList_.push_back(upArcId);
    }

    inline void addUpSuperArcId(const int &upSuperArcId) {
      upSuperArcList_.push_back(upSuperArcId);
    }

    inline int getDownArcId(const int &neighborId) const {
      if((neighborId < 0) || (neighborId >= (int)downArcList_.size()))
        return -1;
      return downArcList_[neighborId];
    }

    inline int getDownSuperArcId(const int &neighborId) const {
      if((neighborId < 0) || (neighborId >= (int)downSuperArcList_.size()))
        return -1;
      return downSuperArcList_[neighborId];
    }

    inline int getUpArcId(const int &neighborId) const {
      if((neighborId < 0) || (neighborId >= (int)upArcList_.size()))
        return -1;
      return upArcList_[neighborId];
    }

    inline int getUpSuperArcId(const int &neighborId) const {
      if((neighborId < 0) || (neighborId >= (int)upSuperArcList_.size()))
        return -1;
      return upSuperArcList_[neighborId];
    }

    inline int getNumberOfDownArcs() const {
      return (int)downArcList_.size();
    }

    inline int getNumberOfDownSuperArcs() const {
      return (int)downSuperArcList_.size();
    }

    inline int getNumberOfUpArcs() const {
      return (int)upArcList_.size();
    }

    inline int getNumberOfUpSuperArcs() const {
      return (int)upSuperArcList_.size();
    }

    inline int getVertexId() const {
      return vertexId_;
    }

    inline int removeDownArcId(const int &arcId) {
      if((arcId < 0) || (arcId >= (int)downArcList_.size()))
        return -1;
      downArcList_.erase(downArcList_.begin() + arcId);
      return 0;
    }

    inline int removeUpArcId(const int &arcId) {
      if((arcId < 0) || (arcId >= (int)upArcList_.size()))
        return -1;
      upArcList_.erase(upArcList_.begin() + arcId);
      return 0;
    }

    inline void setVertexId(const int &vertexId) {
      vertexId_ = vertexId;
    }

  protected:
    friend class SubLevelSetTree;

    int vertexId_{-1};
    bool pruned_{false};
    double layoutX_{}, layoutY_{};
    std::vector<int> downArcList_{}, upArcList_{};
    std::vector<int> downSuperArcList_{}, upSuperArcList_{};
  };

  class ContourTreeSimplificationMetric : virtual public Debug {
  public:
    friend class SubLevelSetTree;

    virtual double
      computeSuperArcMetric(const int &downVertexId,
                            const int &upVertexId,
                            const std::vector<int> &interiorNodeIds)
      = 0;

  protected:
    SubLevelSetTree *tree_{};
  };

  class PersistenceMetric : public ContourTreeSimplificationMetric {
  public:
    double
      computeSuperArcMetric(const int &downVertexId,
                            const int &upVertexId,
                            const std::vector<int> &interiorNodeIds) override;
  };

  class SubLevelSetTree : virtual public Debug {
  public:
    SubLevelSetTree();

    int build();

    // the output list is sorted in ascending (respectively descending)
    // order of function value for the merge tree (respectively for the split
    // tree)
    int buildExtremumList(std::vector<int> &extremumList,
                          const bool &isSubLevelSet = true);

    bool buildPlanarLayout(const double &scaleX, const double &scaleY);

    int buildSaddleList(std::vector<int> &vertexList) const;

    int clearArc(const int &vertexId0, const int &vertexId1);

    int clearRegularNode(const int &vertexId);

    int clearRoot(const int &vertexId);

    int exportPersistenceCurve(const std::string &fileName
                               = "output.plot") const;

    int exportPersistenceDiagram(const std::string &fileName
                                 = "output.plot") const;

    int exportToSvg(const std::string &fileName,
                    const double &scaleX = 1,
                    const double &scaleY = 1);

    int exportToVtk(const std::string &fileName,
                    // fixes a bug in paraview, the voxel size of the cube file
                    // format is not taken into account...
                    const std::vector<float> *origin = NULL,
                    const std::vector<float> *voxelSize = NULL);

    int flush();

    inline const Arc *getArc(const int &arcId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((arcId < 0) || (arcId >= (int)arcList_.size()))
        this->printErr("Out-of-bounds access in getArc: element "
                       + std::to_string(arcId) + " in list of size "
                       + std::to_string(arcList_.size()));
#endif // !TTK_ENABLE_KAMIKAZE
      return &(arcList_[arcId]);
    }

    // this list is sorted only if buildExtremumList has been called before.
    inline const std::vector<int> *getExtremumList() const {
      if(minimumList_)
        return minimumList_;
      return maximumList_;
    }

    inline const Node *getNode(const int &nodeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
        this->printErr("Out-of-bounds access in getNode: element "
                       + std::to_string(nodeId) + " in list of size "
                       + std::to_string(nodeList_.size()));
#endif // !TTK_ENABLE_KAMIKAZE
      return &(nodeList_[nodeId]);
    }

    inline const Node *getNodeDownNeighbor(const Node *n,
                                           const int &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(n == nullptr)
        this->printErr("Nullptr dereference in getNodeDownNeighbor");
#endif // !TTK_ENABLE_KAMIKAZE
      return getNodeDownNeighbor(n - &(nodeList_[0]), neighborId);
    }

    inline const Node *getNodeDownNeighbor(const int &nodeId,
                                           const int &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
        this->printErr("Out-of-bounds access in getNodeDownNeighbor: element "
                       + std::to_string(nodeId) + " in list of size "
                       + std::to_string(nodeList_.size()));
#endif // !TTK_ENABLE_KAMIKAZE
      return &(nodeList_[arcList_[nodeList_[nodeId].getDownArcId(neighborId)]
                           .getDownNodeId()]);
    }

    inline const Node *getNodeUpNeighbor(const Node *n,
                                         const int &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(n == nullptr)
        this->printErr("Nullptr dereference in getNodeUpNeighbor");
#endif // !TTK_ENABLE_KAMIKAZE
      return getNodeUpNeighbor(n - &(nodeList_[0]), neighborId);
    }

    inline const Node *getNodeUpNeighbor(const int &nodeId,
                                         const int &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
        this->printErr("Out-of-bounds access in getNodeUpNeighbor: element "
                       + std::to_string(nodeId) + " in list of size "
                       + std::to_string(nodeList_.size()));
#endif // !TTK_ENABLE_KAMIKAZE
      return &(nodeList_[arcList_[nodeList_[nodeId].getUpArcId(neighborId)]
                           .getUpNodeId()]);
    }

    inline double getNodeScalar(const int &nodeId) const {
      if(!vertexScalars_)
        return -DBL_MAX;
      if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
        return -DBL_MAX;
      return (*vertexScalars_)[nodeList_[nodeId].vertexId_];
    }

    inline int getNumberOfArcs() const {
      return (int)arcList_.size();
    }

    inline int getNumberOfSuperArcs() const {
      return (int)superArcList_.size();
    }

    inline int getNumberOfNodes() const {
      return (int)nodeList_.size();
    }

    int getPersistenceDiagram(
      std::vector<std::pair<double, double>> &diagram,
      std::vector<std::pair<std::pair<int, int>, double>> *pairs
      = nullptr) const;

    // std::vector:
    // - vertex std::pair:
    //   - extremum vertex Id (first.first)
    //   - std::paired saddle Id (first.second)
    // - persistence (second)
    virtual int getPersistencePairs(
      std::vector<std::pair<std::pair<int, int>, double>> &pairs,
      std::vector<std::pair<std::pair<int, int>, double>> *mergePairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *splitPairs
      = nullptr) const;

    int getPersistencePlot(
      std::vector<std::pair<double, int>> &plot,
      std::vector<std::pair<std::pair<int, int>, double>> *persistencePairs
      = nullptr) const;

    inline const SuperArc *getSuperArc(const int &superArcId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((superArcId < 0) || (superArcId >= (int)superArcList_.size()))
        this->printErr("Out-of-bounds access in getSuperArc: element "
                       + std::to_string(superArcId) + " in list of size "
                       + std::to_string(superArcList_.size()));
#endif // !TTK_ENABLE_KAMIKAZE
      return &(superArcList_[superArcId]);
    }

    inline int getVertexScalar(const int &vertexId, double &scalar) {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!vertexScalars_)
        return -1;
      if((vertexId < 0) || (vertexId >= (int)vertexScalars_->size()))
        return -2;
#endif

      scalar = (*vertexScalars_)[vertexId];

      return 0;
    }

    inline const SuperArc *getVertexSuperArc(const int &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        this->printErr("Out-of-bounds access in getVertexSuperArc: element "
                       + std::to_string(vertexId) + " in list of size "
                       + std::to_string(vertexNumber_));
      if(vertex2superArc_[vertexId] == -1)
        this->printErr("Invalid super arc id for vertex "
                       + std::to_string(vertexId));
#endif // !TTK_ENABLE_KAMIKAZE

      return &(superArcList_[vertex2superArc_[vertexId]]);
    }

    inline int getVertexSuperArcId(const int &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        this->printErr("Out-of-bounds access in getVertexSuperArcId: element "
                       + std::to_string(vertexId) + " in list of size "
                       + std::to_string(vertexNumber_));
#endif // !TTK_ENABLE_KAMIKAZE
      return vertex2superArc_[vertexId];
    }

    inline const Node *getVertexNode(const int &vertexId) const {
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return NULL;
      if(vertex2node_[vertexId] != -1)
        return &(nodeList_[vertex2node_[vertexId]]);
      return NULL;
    }

    inline int getVertexNodeId(const int &vertexId) const {
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      return vertex2node_[vertexId];
    }

    bool isJoinTree() const {
      return ((maximumList_ == NULL)
              || ((maximumList_) && (maximumList_->empty())));
    }

    bool isSplitTree() const {
      return ((minimumList_ == NULL)
              || ((minimumList_) && (minimumList_->empty())));
    }

    bool isSosLowerThan(const int &vertexId0, const int &vertexId1) const;

    bool isSosHigherThan(const int &vertexId0, const int &vertexId1) const;

    virtual inline int maintainRegularVertices(const bool &onOff) {
      maintainRegularVertices_ = onOff;
      return 0;
    }

    int moveRegularNode(const Node *n,
                        const Node *oldDown,
                        const Node *oldUp,
                        const Node *newDown,
                        const Node *newUp);

    int print() const;

    inline void setMaximumList(std::vector<int> &maximumList) {
      maximumList_ = &(maximumList);
    }

    inline void setMinimumList(std::vector<int> &minimumList) {
      minimumList_ = &(minimumList);
    }

    inline void setNumberOfVertices(const int &vertexNumber) {
      vertexNumber_ = vertexNumber;
      vertex2superArc_.resize(vertexNumber, -1);
      vertex2superArcNode_.resize(vertexNumber, -1);
      vertex2node_.resize(vertexNumber, -1);
    }

    inline void
      setTriangulation(const AbstractTriangulation *const triangulation) {
      triangulation_ = triangulation;
    }

    inline void
      setVertexPositions(std::vector<std::vector<double>> *vertexPositions) {
      vertexPositions_ = vertexPositions;
    }

    inline void setVertexScalars(const std::vector<real> *const vertexScalars) {
      vertexScalars_ = vertexScalars;
      minScalar_ = 0, maxScalar_ = 0;
      for(int i = 0; i < (int)vertexScalars_->size(); i++) {
        if((!i) || (minScalar_ > (*vertexScalars_)[i])) {
          minScalar_ = (*vertexScalars_)[i];
        }
        if((!i) || (maxScalar_ < (*vertexScalars_)[i])) {
          maxScalar_ = (*vertexScalars_)[i];
        }
      }
    }

    inline void setVertexSoSoffsets(std::vector<int> *vertexSoSoffsets) {
      vertexSoSoffsets_ = vertexSoSoffsets;
    }

    virtual int simplify(const double &simplificationThreshold,
                         ContourTreeSimplificationMetric *metric = NULL);

    int sample(unsigned int samplingLevel = 3);
    int computeBarycenters();

    int getSkeletonScalars(
      const std::vector<double> &scalars,
      std::vector<std::vector<double>> &skeletonScalars) const;
    virtual int computeSkeleton(unsigned int arcResolution = 3);
    virtual int smoothSkeleton(unsigned int skeletonSmoothing);
    virtual int clearSkeleton();

  protected:
    int appendRegularNode(const int &superArcId, const int &nodeId);

    int closeSuperArc(const int &superArcId, const int &nodeId);

    int exportNodeColorToVtk(const int &nodeId, std::ofstream &o);

    int exportNodePosToVtk(const int &nodeId,
                           const int &pointId,
                           std::vector<int> &vertexIds,
                           const std::vector<float> *origin,
                           const std::vector<float> *voxelSize,
                           std::ofstream &o);

    int exportArcPosToVtk(const int &arcId,
                          const int &pointId,
                          std::vector<int> &vertexIds,
                          const std::vector<float> *origin,
                          const std::vector<float> *voxelSize,
                          std::ofstream &o);

    int makeArc(const int &nodeId0, const int &nodeId1);

    int makeNode(const int &vertexId);

    int openSuperArc(const int &nodeId);

    int vertexNumber_{0};
    bool maintainRegularVertices_{true};
    double minScalar_{}, maxScalar_{};
    const std::vector<real> *vertexScalars_{};
    std::vector<int> *vertexSoSoffsets_{};
    const AbstractTriangulation *triangulation_{};
    std::vector<int> *minimumList_{}, *maximumList_{};
    std::vector<Node> nodeList_{}, originalNodeList_{};
    std::vector<Arc> arcList_{};
    std::vector<SuperArc> superArcList_{}, originalSuperArcList_{};
    std::vector<int> vertex2node_{}, vertex2superArc_{}, vertex2superArcNode_{};
    std::vector<std::vector<double>> *vertexPositions_{};
    bool isSkeletonComputed_{false};
  };

  class ContourTree : public SubLevelSetTree {
  public:
    ContourTree();

    int build();

    inline const SubLevelSetTree *getMergeTree() const {
      return &mergeTree_;
    }

    int getPersistencePairs(
      std::vector<std::pair<std::pair<int, int>, double>> &pairs,
      std::vector<std::pair<std::pair<int, int>, double>> *mergePairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *splitPairs
      = nullptr) const;

    int getPersistencePlot(
      std::vector<std::pair<double, int>> &plot,
      std::vector<std::pair<std::pair<int, int>, double>> *mergePairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *splitPairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *pairs
      = nullptr) const;

    int getPersistenceDiagram(
      std::vector<std::pair<double, double>> &diagram,
      std::vector<std::pair<std::pair<int, int>, double>> *mergePairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *splitPairs = nullptr,
      std::vector<std::pair<std::pair<int, int>, double>> *pairs
      = nullptr) const;

    inline const SubLevelSetTree *getSplitTree() const {
      return &splitTree_;
    }

    inline int maintainRegularVertices(const bool &onOff) {
      mergeTree_.maintainRegularVertices(onOff);
      splitTree_.maintainRegularVertices(onOff);
      return 0;
    }

    int
      setVertexNeighbors(const std::vector<std::vector<int>> *vertexNeighbors);

    int setVertexNeighbors(const int &vertexId,
                           const std::vector<int> &neighborList);

    int computeSkeleton(unsigned int arcResolution = 3);
    int smoothSkeleton(unsigned int skeletonSmoothing);
    int clearSkeleton();

    int simplify(const double &simplificationThreshold,
                 ContourTreeSimplificationMetric *metric = NULL);

  protected:
    int combineTrees();

    int finalize();

    int finalizeSuperArc(const int &nodeId, const int &arcId);

    bool isNodeEligible(const Node *n) const;

    SubLevelSetTree mergeTree_{}, splitTree_{};
  };
} // namespace ttk
