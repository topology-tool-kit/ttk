/// \ingroup base
/// \class ttk::CompactTriangulation
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date November 2021.
///
/// \brief CompactTriangulation is a class implemented based on the
/// TopoCluster data structure, which is a localized data structure for
/// simplicial meshes.
/// The key idea of TopoCluster is to subdivide the simplicial mesh into
/// clusters. Then, the connectivity information is computed locally for each
/// cluster and discarded when it is no longer needed. The simplicial mesh needs
/// to be subdivided before using TopoCluster, i.e., there has to be a scalar
/// field named "ttkCompactTriangulationIndex" to denote the cluster index of
/// each vertex. Note Topocluster will reindex the simplices based on the
/// clustering input array.
/// \b Related \b publications \n
/// "TopoCluster: A Localized Data Structure for Topology-based Visualization"
/// Guoxi Liu, Federico Iuricich, Riccardo Fellegara, and Leila De Floriani
/// IEEE Transactions on Visualization and Computer Graphics, 2021.
///
/// \sa ttk::Triangulation
/// \sa ttk::CompactTriangulationPreconditioning

#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CellArray.h>
#include <FlatJaggedArray.h>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <list>

namespace ttk {

  class ImplicitCluster {
  private:
    /* components */
    SimplexId nid;
    std::vector<std::array<SimplexId, 2>> internalEdgeList_;
    std::vector<std::array<SimplexId, 3>> internalTriangleList_;
    boost::unordered_map<std::array<SimplexId, 2>, SimplexId> internalEdgeMap_;
    boost::unordered_map<std::array<SimplexId, 2>, SimplexId> externalEdgeMap_;
    boost::unordered_map<std::array<SimplexId, 3>, SimplexId>
      internalTriangleMap_;
    boost::unordered_map<std::array<SimplexId, 3>, SimplexId>
      externalTriangleMap_;
    /* boundary cells */
    std::vector<bool> boundaryVertices_;
    std::vector<bool> boundaryEdges_;
    std::vector<bool> boundaryTriangles_;
    /* vertex relationships */
    FlatJaggedArray vertexEdges_;
    FlatJaggedArray vertexLinks_;
    FlatJaggedArray vertexNeighbors_;
    FlatJaggedArray vertexStars_;
    FlatJaggedArray vertexTriangles_;
    /* edge relationships */
    // edgeVertex relation can be extracted from internal edge list
    FlatJaggedArray edgeLinks_;
    FlatJaggedArray edgeStars_;
    FlatJaggedArray edgeTriangles_;
    /* triangle relationships */
    // triangleVertex relation can be extracted from internal triangle list
    std::vector<std::array<SimplexId, 3>> triangleEdges_;
    FlatJaggedArray triangleLinks_;
    FlatJaggedArray triangleStars_;
    /* cell relationships */
    std::vector<std::array<SimplexId, 6>> tetraEdges_;
    FlatJaggedArray cellNeighbors_;
    std::vector<std::array<SimplexId, 4>> tetraTriangles_;

  public:
    ImplicitCluster() {
    }
    ImplicitCluster(SimplexId id) : nid(id) {
    }
    ~ImplicitCluster() {
    }

    friend class CompactTriangulation;
  };

  class CompactTriangulation final : public AbstractTriangulation {

    // different id types for compact triangulation
    enum class SIMPLEX_ID { EDGE_ID = 1, TRIANGLE_ID = 2 };

  public:
    CompactTriangulation();
    CompactTriangulation(const CompactTriangulation &rhs);
    CompactTriangulation &operator=(const CompactTriangulation &rhs);
    ~CompactTriangulation();
    CompactTriangulation(CompactTriangulation &&rhs) = default;
    CompactTriangulation &operator=(CompactTriangulation &&rhs) = default;

    /**
     * Set up vertices from the input.
     */
    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const int *indexArray,
                              const bool &doublePrecision = false) {

      if(vertexNumber_)
        clear();

      vertexNumber_ = pointNumber;
      pointSet_ = pointSet;
      vertexIndices_ = indexArray;
      doublePrecision_ = doublePrecision;

      return 0;
    }

    /**
     * Set up cells from the input.
     */

#ifdef TTK_CELL_ARRAY_NEW
    // Layout with connectivity + offset array (new)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *connectivity,
                             const LongSimplexId *offset) {

      // Cell Check
      {
        if(cellNumber > 0) {
          const auto &cellDimension = offset[1] - offset[0] - 1;

          if(cellDimension < 0 || cellDimension > 3) {
            this->printErr("Unable to create triangulation for cells of "
                           "dimension 4 or higher ("
                           + std::to_string(cellDimension) + ").");
            return -1;
          }

          bool error = false;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
          for(SimplexId i = 0; i < cellNumber; i++) {
            if(offset[i + 1] - offset[i] - 1 != cellDimension) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
              error = true;
            }
          }

          if(error) {
            this->printErr("Unable to create triangulation for "
                           "inhomogeneous\ncell dimensions.");
            return -2;
          }
        }
      }

      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;
      std::vector<SimplexId> vertexMap(vertexNumber_);
      reorderVertices(vertexMap);
      reorderCells(vertexMap, cellNumber, connectivity, offset);
      cellArray_
        = std::make_shared<CellArray>(connectivity, offset, cellNumber);

      // ASSUME Regular Mesh Here to compute dimension!
      if(cellNumber) {
        if(cellArray_->getCellVertexNumber(0) == 3) {
          maxCellDim_ = 2;
          triangleIntervals_ = cellIntervals_;
        } else {
          maxCellDim_ = 3;
        }
      }

      return 0;
    }
#else
    // Flat layout with a single array (legacy & default one)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {
      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;

      if(cellNumber) {
        // assume regular mesh here to compute dimension
        maxCellDim_ = cellArray[0] - 1;
        std::vector<SimplexId> vertexMap(vertexNumber_);
        reorderVertices(vertexMap);
        reorderCells(vertexMap, cellArray);
        cellArray_ = std::make_shared<CellArray>(
          cellArray, cellNumber, cellArray[0] - 1);
      }

      return 0;
    }
#endif

    /**
     * Reorder the input vertices.
     */
    int reorderVertices(std::vector<SimplexId> &vertexMap);

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const SimplexId &cellNumber,
                     const LongSimplexId *connectivity,
                     const LongSimplexId *offset);

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const LongSimplexId *cellArray);

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_)) {
        edgeId = -1;
        return 0;
      }
      if(localEdgeId < 0) {
        edgeId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->tetraEdges_.empty()) {
        getClusterTetraEdges(exnode);
      }

      if(localEdgeId >= (int)(exnode->tetraEdges_)[localCellId].size()) {
        edgeId = -2;
      } else {
        edgeId = (exnode->tetraEdges_)[localCellId][localEdgeId];
      }
      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      (void)cellId;
      return (maxCellDim_ + 1) * maxCellDim_ / 2;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() override {
      if(cellEdgeVector_.empty()) {
        cellEdgeVector_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          ImplicitCluster *exnode = searchCache(nid);
          if(exnode->tetraEdges_.empty()) {
            getClusterTetraEdges(exnode);
          }
          for(size_t i = 0; i < exnode->tetraEdges_.size(); i++) {
            cellEdgeVector_.push_back({exnode->tetraEdges_.at(i).begin(),
                                       exnode->tetraEdges_.at(i).end()});
          }
        }
      }
      return &cellEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_)) {
        neighborId = -1;
        return 0;
      }
      if(localNeighborId < 0) {
        neighborId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->cellNeighbors_.empty()) {
        getClusterCellNeighbors(exnode);
      }

      if(localNeighborId >= exnode->cellNeighbors_.size(localCellId)) {
        neighborId = -2;
      } else {
        neighborId = exnode->cellNeighbors_.get(localCellId, localNeighborId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->cellNeighbors_.empty()) {
        getClusterCellNeighbors(exnode);
      }
      return exnode->cellNeighbors_.size(localCellId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborList_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localCellNeighbors;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->cellNeighbors_.empty()) {
          getClusterCellNeighbors(exnode);
        }
        exnode->cellNeighbors_.copyTo(localCellNeighbors);
        cellNeighborList_.insert(cellNeighborList_.end(),
                                 localCellNeighbors.begin(),
                                 localCellNeighbors.end());
      }
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_)) {
        triangleId = -1;
        return 0;
      }
      if((localTriangleId < 0)
         || (localTriangleId >= getCellTriangleNumber(cellId))) {
        triangleId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->tetraTriangles_.empty()) {
        getClusterCellTriangles(exnode);
      }
      triangleId = (exnode->tetraTriangles_)[localCellId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      (void)cellId;
      return (maxCellDim_ + 1) * maxCellDim_ * (maxCellDim_ - 1) / 6;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override {
      if(cellTriangleVector_.empty()) {
        cellTriangleVector_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          ImplicitCluster *exnode = searchCache(nid);
          if(exnode->tetraTriangles_.empty()) {
            getClusterCellTriangles(exnode);
          }
          for(size_t i = 0; i < exnode->tetraTriangles_.size(); i++) {
            cellTriangleVector_.push_back(
              {exnode->tetraTriangles_.at(i).begin(),
               exnode->tetraTriangles_.at(i).end()});
          }
        }
      }
      return &cellTriangleVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_)) {
        vertexId = -1;
        return 0;
      }
      if((localVertexId < 0)
         || (localVertexId >= cellArray_->getCellVertexNumber(cellId))) {
        vertexId = -2;
        return 0;
      }
#endif

      vertexId = cellArray_->getCellVertex(cellId, localVertexId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      return cellArray_->getCellVertexNumber(cellId);
    }

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return maxCellDim_;
    }

    inline const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      edgeList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->internalEdgeList_.empty())
          buildInternalEdgeMap(exnode, true, false);
        edgeList_.insert(edgeList_.end(), exnode->internalEdgeList_.begin(),
                         exnode->internalEdgeList_.end());
      }
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back())) {
        linkId = -1;
        return 0;
      }
      if(localLinkId < 0) {
        linkId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeLinks_.empty()) {
        getClusterEdgeLinks(exnode);
      }

      if(localLinkId >= exnode->edgeLinks_.size(localEdgeId)) {
        linkId = -2;
      } else {
        linkId = exnode->edgeLinks_.get(localEdgeId, localLinkId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeLinks_.empty()) {
        getClusterEdgeLinks(exnode);
      }
      return exnode->edgeLinks_.size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeLinks;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeLinks_.empty()) {
          getClusterEdgeLinks(exnode);
        }
        exnode->edgeLinks_.copyTo(localEdgeLinks);
        edgeLinkList_.insert(
          edgeLinkList_.end(), localEdgeLinks.begin(), localEdgeLinks.end());
      }
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back())) {
        starId = -1;
        return 0;
      }
      if(localStarId < 0) {
        starId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeStars_.empty()) {
        getClusterEdgeStars(exnode);
      }

      if(localStarId >= exnode->edgeStars_.size(localEdgeId)) {
        starId = -2;
      } else {
        starId = exnode->edgeStars_.get(localEdgeId, localStarId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
#endif
      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      ImplicitCluster *exnode = searchCache(nid);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      if(exnode->edgeStars_.empty()) {
        getClusterEdgeStars(exnode);
      }
      return exnode->edgeStars_.size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeStars;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeStars_.empty()) {
          getClusterEdgeStars(exnode);
        }
        exnode->edgeStars_.copyTo(localEdgeStars);
        edgeStarList_.insert(
          edgeStarList_.end(), localEdgeStars.begin(), localEdgeStars.end());
      }
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back())) {
        triangleId = -1;
        return 0;
      }
      if(localTriangleId < 0) {
        triangleId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeTriangles_.empty()) {
        getClusterEdgeTriangles(exnode);
      }

      if(localTriangleId >= exnode->edgeTriangles_.size(localEdgeId)) {
        triangleId = -2;
      } else {
        triangleId = exnode->edgeTriangles_.get(localEdgeId, localTriangleId);
      }
      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeTriangles_.empty()) {
        getClusterEdgeTriangles(exnode);
      }
      return exnode->edgeTriangles_.size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeTriangles;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeTriangles_.empty()) {
          getClusterEdgeTriangles(exnode);
        }
        exnode->edgeTriangles_.copyTo(localEdgeTriangles);
        edgeTriangleList_.insert(edgeTriangleList_.end(),
                                 localEdgeTriangles.begin(),
                                 localEdgeTriangles.end());
      }
      return &edgeTriangleList_;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back())) {
        vertexId = -1;
        return 0;
      }
      if((localVertexId != 0) && (localVertexId != 1)) {
        vertexId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->internalEdgeList_.empty()) {
        buildInternalEdgeMap(exnode, true, false);
      }

      if(localVertexId) {
        vertexId = exnode->internalEdgeList_.at(localEdgeId)[1];
      } else {
        vertexId = exnode->internalEdgeList_.at(localEdgeId)[0];
      }
      return 0;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdgesInternal() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!edgeIntervals_.size())
        return -1;
#endif

      return edgeIntervals_.back() + 1;
    }

    inline SimplexId getNumberOfTrianglesInternal() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!triangleIntervals_.size())
        return -1;
#endif

      return triangleIntervals_.back() + 1;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    inline const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override {
      // if it is a triangle mesh
      if(getDimensionality() == 2) {
        triangleList_.resize(cellNumber_, std::array<SimplexId, 3>());
        for(SimplexId cid = 0; cid < cellNumber_; cid++) {
          triangleList_[cid][0] = cellArray_->getCellVertex(cid, 0);
          triangleList_[cid][1] = cellArray_->getCellVertex(cid, 1);
          triangleList_[cid][2] = cellArray_->getCellVertex(cid, 2);
        }
      } else {
        triangleList_.reserve(triangleIntervals_.back() + 1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          ImplicitCluster *exnode = searchCache(nid);
          if(exnode->internalTriangleList_.empty()) {
            buildInternalTriangleMap(exnode, true, false);
          }
          triangleList_.insert(triangleList_.end(),
                               exnode->internalTriangleList_.begin(),
                               exnode->internalTriangleList_.end());
        }
      }
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 || triangleId > triangleIntervals_.back()) {
        edgeId = -1;
        return 0;
      }
      if((localEdgeId < 0) || (localEdgeId > 2)) {
        edgeId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleEdges_.empty()) {
        getClusterTriangleEdges(exnode);
      }
      edgeId = (exnode->triangleEdges_)[localTriangleId][localEdgeId];
      return 0;
    }

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      (void)triangleId;
      return 3;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override {
      if(triangleEdgeVector_.empty()) {
        triangleEdgeVector_.reserve(triangleIntervals_.size() + 1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          ImplicitCluster *exnode = searchCache(nid);
          if(exnode->triangleEdges_.empty()) {
            getClusterTriangleEdges(exnode);
          }
          for(size_t i = 0; i < exnode->triangleEdges_.size(); i++) {
            triangleEdgeVector_.push_back({exnode->triangleEdges_.at(i).begin(),
                                           exnode->triangleEdges_.at(i).end()});
          }
        }
      }
      return &triangleEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back())) {
        linkId = -1;
        return 0;
      }
      if(localLinkId < 0) {
        linkId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleLinks_.empty()) {
        getClusterTriangleLinks(exnode);
      }

      if(localLinkId >= exnode->triangleLinks_.size(localTriangleId)) {
        linkId = -2;
      } else {
        linkId = exnode->triangleLinks_.get(localTriangleId, localLinkId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleLinks_.empty()) {
        getClusterTriangleLinks(exnode);
      }
      return exnode->triangleLinks_.size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleLinks;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->triangleLinks_.empty()) {
          getClusterTriangleLinks(exnode);
        }
        exnode->triangleLinks_.copyTo(localTriangleLinks);
        triangleLinkList_.insert(triangleLinkList_.end(),
                                 localTriangleLinks.begin(),
                                 localTriangleLinks.end());
      }
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || !triangleIntervals_.size()
         || (triangleId > triangleIntervals_.back())) {
        starId = -1;
        return 0;
      }
      if(localStarId < 0) {
        starId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleStars_.empty()) {
        getClusterTriangleStars(exnode);
      }

      if(localStarId >= exnode->triangleStars_.size(localTriangleId)) {
        starId = -2;
      } else {
        starId = exnode->triangleStars_.get(localTriangleId, localStarId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || !triangleIntervals_.size()
         || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleStars_.empty()) {
        getClusterTriangleStars(exnode);
      }
      return exnode->triangleStars_.size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleStars;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->triangleStars_.empty()) {
          getClusterTriangleStars(exnode);
        }
        triangleStarList_.insert(triangleStarList_.end(),
                                 localTriangleStars.begin(),
                                 localTriangleStars.end());
      }
      return &triangleStarList_;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back())) {
        vertexId = -1;
        return 0;
      }
      if((localVertexId < 0) || (localVertexId > 2)) {
        vertexId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->internalTriangleList_.empty()) {
        buildInternalTriangleMap(exnode, true, false);
      }
      vertexId
        = exnode->internalTriangleList_.at(localTriangleId)[localVertexId];
      return 0;
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_)) {
        edgeId = -1;
        return 0;
      }
      if(localEdgeId < 0) {
        edgeId = -1;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexEdges_.empty()) {
        getClusterVertexEdges(exnode);
      }
      if(localEdgeId >= exnode->vertexEdges_.size(localVertexId)) {
        edgeId = -2;
      } else {
        edgeId = exnode->vertexEdges_.get(localVertexId, localEdgeId);
      }
      return 0;
    }

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexEdges_.empty()) {
        getClusterVertexEdges(exnode);
      }
      return exnode->vertexEdges_.size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      vertexEdgeList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexEdges;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexEdges_.empty()) {
          getClusterVertexEdges(exnode);
        }
        exnode->vertexEdges_.copyTo(localVertexEdges);
        vertexEdgeList_.insert(vertexEdgeList_.end(), localVertexEdges.begin(),
                               localVertexEdges.end());
      }

      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId > vertexIntervals_.back())) {
        linkId = -1;
        return 0;
      }
      if(localLinkId < 0) {
        linkId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexLinks_.empty()) {
        getClusterVertexLinks(exnode);
      }
      if(localLinkId >= exnode->vertexLinks_.size(localVertexId)) {
        linkId = -2;
      } else {
        linkId = exnode->vertexLinks_.get(localVertexId, localLinkId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId > vertexIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexLinks_.empty()) {
        getClusterVertexLinks(exnode);
      }
      return exnode->vertexLinks_.size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkList_.reserve(vertexIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexLinks;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexLinks_.empty()) {
          getClusterVertexLinks(exnode);
        }
        exnode->vertexLinks_.copyTo(localVertexLinks);
        vertexLinkList_.insert(vertexLinkList_.end(), localVertexLinks.begin(),
                               localVertexLinks.end());
      }
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_)) {
        neighborId = -1;
        return 0;
      }
      if(localNeighborId < 0) {
        neighborId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexNeighbors_.empty()) {
        getClusterVertexNeighbors(exnode);
      }
      if(localNeighborId >= exnode->vertexNeighbors_.size(localVertexId)) {
        neighborId = -2;
      } else {
        neighborId
          = exnode->vertexNeighbors_.get(localVertexId, localNeighborId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexNeighbors_.empty()) {
        getClusterVertexNeighbors(exnode);
      }
      return exnode->vertexNeighbors_.size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexNeighbors;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexNeighbors_.empty()) {
          getClusterVertexNeighbors(exnode);
        }
        exnode->vertexNeighbors_.copyTo(localVertexNeighbors);
        vertexNeighborList_.insert(vertexNeighborList_.end(),
                                   localVertexNeighbors.begin(),
                                   localVertexNeighbors.end());
      }
      return &vertexNeighborList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
      const SimplexId &vertexId, float &x, float &y, float &z) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      if(doublePrecision_) {
        x = ((double *)pointSet_)[3 * vertexId];
        y = ((double *)pointSet_)[3 * vertexId + 1];
        z = ((double *)pointSet_)[3 * vertexId + 2];
      } else {
        x = ((float *)pointSet_)[3 * vertexId];
        y = ((float *)pointSet_)[3 * vertexId + 1];
        z = ((float *)pointSet_)[3 * vertexId + 2];
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_)) {
        starId = -1;
        return 0;
      }
      if(localStarId < 0) {
        starId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexStars_.empty()) {
        getClusterVertexStars(exnode);
      }
      if(localStarId >= exnode->vertexStars_.size(localVertexId)) {
        starId = -2;
      } else {
        starId = exnode->vertexStars_.get(localVertexId, localStarId);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexStars_.empty()) {
        getClusterVertexStars(exnode);
      }
      return exnode->vertexStars_.size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexStars;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexStars_.empty()) {
          getClusterVertexStars(exnode);
        }
        exnode->vertexStars_.copyTo(localVertexStars);
        vertexStarList_.insert(vertexStarList_.end(), localVertexStars.begin(),
                               localVertexStars.end());
      }
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_)) {
        triangleId = -1;
        return 0;
      }
      if(localTriangleId < 0) {
        triangleId = -2;
        return 0;
      }
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexTriangles_.empty()) {
        getClusterVertexTriangles(exnode);
      }
      if(localTriangleId >= exnode->vertexTriangles_.size(localVertexId)) {
        triangleId = -2;
      } else {
        triangleId
          = exnode->vertexTriangles_.get(localVertexId, localTriangleId);
      }
      return 0;
    }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexTriangles_.empty()) {
        getClusterVertexTriangles(exnode);
      }
      return exnode->vertexTriangles_.size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexTriangles;
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexTriangles_.empty()) {
          getClusterVertexTriangles(exnode);
        }
        exnode->vertexTriangles_.copyTo(localVertexTriangles);
        vertexTriangleList_.insert(vertexTriangleList_.end(),
                                   localVertexTriangles.begin(),
                                   localVertexTriangles.end());
      }
      return &vertexTriangleList_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return false;
#endif
      SimplexId nid = findNodeIndex(edgeId, SIMPLEX_ID::EDGE_ID);
      SimplexId localedgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      getBoundaryCells(exnode, 1);
      return (exnode->boundaryEdges_)[localedgeId];
    }

    bool isEmpty() const override {
      return !vertexNumber_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override {
      if(getDimensionality() == 2)
        return false;

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return false;
#endif
      SimplexId nid = findNodeIndex(triangleId, SIMPLEX_ID::TRIANGLE_ID);
      SimplexId localtriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      getBoundaryCells(exnode);
      return (exnode->boundaryTriangles_)[localtriangleId];
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return false;
#endif
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      getBoundaryCells(exnode, 0);
      return (exnode->boundaryVertices_)[localVertexId];
    }

    inline int preconditionBoundaryEdgesInternal() override {
      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {
      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() override {
      if(getDimensionality() == 2) {
        preconditionEdges();
      } else if(getDimensionality() == 3) {
        preconditionTriangles();
      }
      return 0;
    }

    inline int preconditionCellEdgesInternal() override {
      return 0;
    }

    inline int preconditionCellNeighborsInternal() override {
      return 0;
    }

    inline int preconditionCellTrianglesInternal() override {
      return 0;
    }

    inline int preconditionEdgesInternal() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexNumber_ <= 0)
        return -1;
      if(cellNumber_ <= 0)
        return -2;
      if(!cellArray_)
        return -3;
      if(nodeNumber_ <= 0)
        return -4;
#endif

      if(edgeIntervals_.empty()) {
        Timer t;
        edgeIntervals_.resize(nodeNumber_ + 1);
        edgeIntervals_[0] = -1;
        std::vector<SimplexId> edgeCount(nodeNumber_ + 1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeCount[nid] = countInternalEdges(nid);
        }

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeIntervals_[nid] = edgeIntervals_[nid - 1] + edgeCount[nid];
        }

        this->printMsg("Edges preconditioned in "
                       + std::to_string(t.getElapsedTime()) + " s.");
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() override {
      return 0;
    }

    inline int preconditionEdgeStarsInternal() override {
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() override {
      return 0;
    }

    inline int preconditionTrianglesInternal() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexNumber_ <= 0)
        return -1;
      if(cellNumber_ <= 0)
        return -2;
      if(!cellArray_)
        return -3;
      if(nodeNumber_ <= 0)
        return -4;
#endif

      // build triangle interval list
      if(triangleIntervals_.empty()) {
        Timer t;
        triangleIntervals_.resize(nodeNumber_ + 1);
        triangleIntervals_[0] = -1;
        std::vector<SimplexId> triangleCount(nodeNumber_ + 1);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          triangleCount[nid] = countInternalTriangles(nid);
        }

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          triangleIntervals_[nid]
            = triangleIntervals_[nid - 1] + triangleCount[nid];
        }

        this->printMsg("Triangles preconditioned in "
                       + std::to_string(t.getElapsedTime()) + " s.");
      }

      return 0;
    }

    inline int preconditionTriangleEdgesInternal() override {
      return 0;
    }

    inline int preconditionTriangleLinksInternal() override {
      return 0;
    }

    inline int preconditionTriangleStarsInternal() override {
      return 0;
    }

    inline int preconditionVertexEdgesInternal() override {
      return 0;
    }

    inline int preconditionVertexLinksInternal() override {
      if(getDimensionality() == 2) {
        preconditionEdges();
      } else if(getDimensionality() == 3) {
        preconditionTriangles();
      } else {
        this->printErr("Unsupported dimension for vertex link precondtion");
        return -1;
      }
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() override {
      return 0;
    }

    inline int preconditionVertexStarsInternal() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!cellArray_)
        return -1;
#endif

      return 0;
    }

    inline int preconditionVertexTrianglesInternal() override {
      return 0;
    }

    /**
     * Initialize the cache with the ratio.
     */
    inline void initCache(const float ratio = 0.2) {
      cacheSize_ = nodeNumber_ * ratio + 1;
      for(int i = 0; i < threadNumber_; i++) {
        caches_[i].clear();
        cacheMaps_[i].clear();
      }
      this->printMsg("Initializing cache: " + std::to_string(cacheSize_));
    }

    /**
     * Get the size of the cache.
     */
    inline size_t getCacheSize() {
      return cacheSize_;
    }

    /**
     * Reset the cache size to better fit in the parallel or sequential
     * algorithms.
     */
    inline int resetCache(int option) {
      for(int i = 0; i < threadNumber_; i++) {
        caches_[i].clear();
        cacheMaps_[i].clear();
      }
      if(option) {
        caches_.resize(1);
        cacheMaps_.resize(1);
        cacheSize_ *= threadNumber_;
      } else {
        caches_.resize(threadNumber_);
        cacheMaps_.resize(threadNumber_);
        cacheSize_ /= threadNumber_;
      }
      return 0;
    }

  protected:
    inline int clear() {
      vertexIntervals_.clear();
      edgeIntervals_.clear();
      triangleIntervals_.clear();
      cellIntervals_.clear();
      externalCells_.clear();
      cacheMaps_.clear();
      cacheMaps_.clear();
      return AbstractTriangulation::clear();
    }

    /**
     * Find the corresponding node index given the id.
     */
    inline SimplexId findNodeIndex(SimplexId id, SIMPLEX_ID idType) const {
      const std::vector<SimplexId> *intervals = nullptr;
      // determine which vector to search
      if(idType == SIMPLEX_ID::EDGE_ID) {
        intervals = &edgeIntervals_;
      } else if(idType == SIMPLEX_ID::TRIANGLE_ID) {
        intervals = &triangleIntervals_;
      } else {
        return -1;
      }

      std::vector<SimplexId>::const_iterator low
        = lower_bound(intervals->begin(), intervals->end(), id);
      return (low - intervals->begin());
    }

    /**
     * Search the node in the cache.
     */
    inline ImplicitCluster *searchCache(const SimplexId &nodeId,
                                        const SimplexId reservedId = 0) const {
      ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif

      // cannot find the expanded node in the cache
      if(cacheMaps_[threadId].find(nodeId) == cacheMaps_[threadId].end()) {
        if(caches_[threadId].size() >= cacheSize_) {
          if(caches_[threadId].back().nid == reservedId) {
            return nullptr;
          }
          cacheMaps_[threadId].erase(caches_[threadId].back().nid);
          caches_[threadId].pop_back();
        }
        caches_[threadId].push_front(ImplicitCluster(nodeId));
        cacheMaps_[threadId][nodeId] = caches_[threadId].begin();
      }
      return &(*cacheMaps_[threadId][nodeId]);
    }

    /**
     * Build the internal edge list in the node.
     */
    int buildInternalEdgeMap(ImplicitCluster *const nodePtr,
                             bool computeInternalEdgeList,
                             bool computeInternalEdgeMap) const;
    /**
     * Build the external edge list in the node.
     */
    int buildExternalEdgeMap(ImplicitCluster *const nodePtr) const;

    /**
     * Build the internal triangle list in the node.
     */
    int buildInternalTriangleMap(ImplicitCluster *const nodePtr,
                                 bool computeInternalTriangleList,
                                 bool computeInternalTriangleMap) const;

    /**
     * Build the external triangle list in the node.
     */
    int buildExternalTriangleMap(ImplicitCluster *const nodePtr) const;

    /**
     * Get the number of internal edges in the node.
     */
    SimplexId countInternalEdges(SimplexId nodeId) const;

    /**
     * Get the number of internal triangles in the node.
     */
    int countInternalTriangles(SimplexId nodeId) const;

    /**
     * Get the cell neighbors for all cells in a given cluster.
     */
    int getClusterCellNeighbors(ImplicitCluster *const nodePtr) const;

    /**
     * Get the cell triangles for all cells in a given cluster.
     */
    int getClusterCellTriangles(ImplicitCluster *const nodePtr) const;

    /**
     * Get the edge links for all the edges in a given cluster.
     */
    int getClusterEdgeLinks(ImplicitCluster *const nodePtr) const;

    /**
     * Get the edge stars for all the edges in a given cluster.
     */
    int getClusterEdgeStars(ImplicitCluster *const nodePtr) const;

    /**
     * Get the edge triangles for all the edges in a given cluster.
     */
    int getClusterEdgeTriangles(ImplicitCluster *const nodePtr) const;

    /**
     * Get the edges for all tetrahedra in a given cluster.
     * Check if the tetraEdges_ is NULL before calling the function.
     */
    int getClusterTetraEdges(ImplicitCluster *const nodePtr) const;

    /**
     * Get the triangle edges for all the triangles in a given cluster.
     */
    int getClusterTriangleEdges(ImplicitCluster *const nodePtr) const;

    /**
     * Get the triangle links for all the triangles in a given cluster.
     */
    int getClusterTriangleLinks(ImplicitCluster *const nodePtr) const;

    /**
     * Get the triangle stars for all the triangles in a given cluster.
     */
    int getClusterTriangleStars(ImplicitCluster *const nodePtr) const;

    /**
     * Get the vertex edges for all the vertices in a given cluster.
     */
    int getClusterVertexEdges(ImplicitCluster *const nodePtr) const;

    /**
     * Get the vertex links for all the vertices in a given cluster.
     */
    int getClusterVertexLinks(ImplicitCluster *const nodePtr) const;

    /**
     * Get the vertex neighbors for all the vertices in a given cluster.
     */
    int getClusterVertexNeighbors(ImplicitCluster *const nodePtr) const;

    /**
     * Get the vertex stars for all the vertices in a given cluster.
     * The function is similar as getVertexCells().
     */
    int getClusterVertexStars(ImplicitCluster *const nodePtr) const;

    /**
     * Get the vertex triangles for all the vertices in a given cluster.
     */
    int getClusterVertexTriangles(ImplicitCluster *const nodePtr) const;

    /**
     * Get the boundary cells in a given cluster.
     */
    int getBoundaryCells(ImplicitCluster *const nodePtr,
                         const SimplexId dim = 2) const;

    /**
     * Protected class variables.
     */
    bool doublePrecision_;
    int maxCellDim_;
    SimplexId cellNumber_, vertexNumber_, nodeNumber_;
    const void *pointSet_;
    const int *vertexIndices_;
    std::vector<SimplexId> vertexIntervals_;
    std::vector<SimplexId> edgeIntervals_;
    std::vector<SimplexId> triangleIntervals_;
    std::vector<SimplexId> cellIntervals_;
    std::shared_ptr<CellArray> cellArray_;
    std::vector<std::vector<SimplexId>> externalCells_;

    // Cache system
    size_t cacheSize_;
    mutable std::vector<std::list<ImplicitCluster>> caches_;
    mutable std::vector<
      boost::unordered_map<SimplexId, std::list<ImplicitCluster>::iterator>>
      cacheMaps_;
  };
} // namespace ttk
