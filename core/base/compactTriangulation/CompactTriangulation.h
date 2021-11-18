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

#define EDGE_ID 1
#define TRIANGLE_ID 2

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

  public:
    CompactTriangulation();

    ~CompactTriangulation();

    /**
     * Set up vertices from the input.
     */
    int setInputPoints(const SimplexId &pointNumber,
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
    int reorderVertices(std::vector<SimplexId> &vertexMap) {
      // get the number of nodes (the max value in the array)
      for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
        if(vertexIndices_[vid] > nodeNumber_) {
          nodeNumber_ = vertexIndices_[vid];
        }
      }
      nodeNumber_++; // since the index starts from 0
      std::vector<std::vector<SimplexId>> nodeVertices(nodeNumber_);
      for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
        nodeVertices[vertexIndices_[vid]].push_back(vid);
      }

      // update the vertex intervals
      vertexIntervals_.resize(nodeNumber_ + 1);
      vertexIntervals_[0] = -1;
      SimplexId vertexCount = 0;
      for(SimplexId nid = 0; nid < nodeNumber_; nid++) {
        for(SimplexId vid : nodeVertices[nid]) {
          vertexMap[vid] = vertexCount++;
        }
        vertexIntervals_[nid + 1] = vertexCount - 1;
      }

      // rearange the vertex coordinate values
      if(doublePrecision_) {
        std::vector<double> newPointSet(3 * vertexNumber_);
        for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
          for(int j = 0; j < 3; j++) {
            newPointSet[3 * vertexMap[vid] + j]
              = ((double *)pointSet_)[3 * vid + j];
          }
        }
        for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
          for(int j = 0; j < 3; j++) {
            ((double *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
          }
        }
      } else {
        std::vector<float> newPointSet(3 * vertexNumber_);
        for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
          for(int j = 0; j < 3; j++) {
            newPointSet[3 * vertexMap[vid] + j]
              = ((float *)pointSet_)[3 * vid + j];
          }
        }
        for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
          for(int j = 0; j < 3; j++) {
            ((float *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
          }
        }
      }

      // change the vertex indices
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        for(SimplexId vid = vertexIntervals_[nid - 1] + 1;
            vid <= vertexIntervals_[nid]; vid++) {
          ((int *)vertexIndices_)[vid] = nid;
        }
      }

      return 0;
    }

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const SimplexId &cellNumber,
                     const LongSimplexId *connectivity,
                     const LongSimplexId *offset) {
      // change the indices in cell array
      SimplexId cellCount = 0, verticesPerCell = offset[1] - offset[0];
      std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
      LongSimplexId *cellArr = ((LongSimplexId *)connectivity);

      for(SimplexId cid = 0; cid < cellNumber; cid++) {
        SimplexId cellId = offset[cid];
        for(int j = 0; j < verticesPerCell; j++) {
          cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
        }
        std::sort(cellArr + cellId, cellArr + cellId + verticesPerCell);
        nodeCells[vertexIndices_[cellArr[cellId]]].push_back(cid);
      }

      // rearange the cell array
      cellIntervals_.resize(nodeNumber_ + 1);
      externalCells_.resize(nodeNumber_ + 1);
      cellIntervals_[0] = -1;
      std::vector<LongSimplexId> newCellArray(verticesPerCell * cellNumber);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        for(SimplexId cid : nodeCells[nid]) {
          SimplexId cellId = verticesPerCell * cid;
          SimplexId newCellId = verticesPerCell * cellCount;
          for(int j = 0; j < verticesPerCell; j++) {
            newCellArray[newCellId + j] = connectivity[cellId + j];
            if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
              SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
              if(externalCells_[nodeNum].empty()
                 || externalCells_[nodeNum].back() != cid) {
                externalCells_[nodeNum].push_back(cid);
              }
            }
          }
          cellCount++;
        }
        cellIntervals_[nid] = cellCount - 1;
      }

      // copy the new cell array back to original one
      for(SimplexId i = 0; i < verticesPerCell * cellNumber; i++) {
        ((LongSimplexId *)connectivity)[i] = newCellArray[i];
      }

      return 0;
    }

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const LongSimplexId *cellArray) {
      // change the indices in cell array
      SimplexId cellCount = 0, verticesPerCell = cellArray[0];
      std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
      LongSimplexId *cellArr = ((LongSimplexId *)cellArray);

      for(SimplexId cid = 0; cid < cellNumber_; cid++) {
        SimplexId cellId = (verticesPerCell + 1) * cid;
        for(int j = 1; j <= verticesPerCell; j++) {
          cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
        }
        std::sort(cellArr + cellId + 1, cellArr + cellId + 1 + verticesPerCell);
        nodeCells[vertexIndices_[cellArr[cellId + 1]]].push_back(cid);
      }

      // rearange the cell array
      cellIntervals_.resize(nodeNumber_ + 1);
      externalCells_.resize(nodeNumber_ + 1);
      cellIntervals_[0] = -1;
      std::vector<LongSimplexId> newCellArray((verticesPerCell + 1)
                                              * cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        for(SimplexId cid : nodeCells[nid]) {
          SimplexId cellId = (verticesPerCell + 1) * cid;
          SimplexId newCellId = (verticesPerCell + 1) * cellCount;
          newCellArray[newCellId] = verticesPerCell;
          for(int j = 1; j <= verticesPerCell; j++) {
            newCellArray[newCellId + j] = cellArray[cellId + j];
            if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
              SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
              if(externalCells_[nodeNum].empty()
                 || externalCells_[nodeNum].back() != cid) {
                externalCells_[nodeNum].push_back(cid);
              }
            }
          }
          cellCount++;
        }
        cellIntervals_[nid] = cellCount - 1;
      }

      // copy the new cell array back to original one
      for(SimplexId i = 0; i < (verticesPerCell + 1) * cellNumber_; i++) {
        ((LongSimplexId *)cellArray)[i] = newCellArray[i];
      }

      return 0;
    }

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
        getClusterCellEdges(exnode);
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
            getClusterCellEdges(exnode);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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
      if((triangleId < 0) || (triangleId > triangleIntervals_.back())) {
        edgeId = -1;
        return 0;
      }
      if((localEdgeId < 0) || (localEdgeId > 2)) {
        edgeId = -2;
        return 0;
      }
#endif

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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
      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
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
      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
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
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionEdgesInternal();
        hasPreconditionedBoundaryEdges_ = true;
      } else {
        this->printErr("Unsupported dimension for boundary edge precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionTrianglesInternal();
        hasPreconditionedBoundaryTriangles_ = true;
      } else {
        this->printErr(
          "Unsupported dimension for boundary triangle precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() override {
      if(getDimensionality() == 2) {
        preconditionEdgesInternal();
      } else if(getDimensionality() == 3) {
        preconditionTrianglesInternal();
      }
      hasPreconditionedBoundaryVertices_ = true;
      return 0;
    }

    inline int preconditionCellEdgesInternal() override {
      hasPreconditionedCellEdges_ = true;
      return 0;
    }

    inline int preconditionCellNeighborsInternal() override {
      hasPreconditionedCellNeighbors_ = true;
      return 0;
    }

    inline int preconditionCellTrianglesInternal() override {
      hasPreconditionedCellTriangles_ = true;
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
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionEdges();
        hasPreconditionedEdgeLinks_ = true;
      } else {
        this->printErr("Unsupported dimension for edge link precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionEdgeStarsInternal() override {
      preconditionEdges();
      hasPreconditionedEdgeStars_ = true;
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() override {
      preconditionEdges();
      preconditionTriangles();
      hasPreconditionedEdgeTriangles_ = true;
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

        hasPreconditionedTriangles_ = true;
        this->printMsg("Triangles preconditioned in "
                       + std::to_string(t.getElapsedTime()) + " s.");
      }

      return 0;
    }

    inline int preconditionTriangleEdgesInternal() override {
      preconditionEdges();
      preconditionTriangles();
      hasPreconditionedTriangleEdges_ = true;
      return 0;
    }

    inline int preconditionTriangleLinksInternal() override {
      return 0;
    }

    inline int preconditionTriangleStarsInternal() override {
      preconditionTrianglesInternal();
      return 0;
    }

    inline int preconditionVertexEdgesInternal() override {
      return 0;
    }

    inline int preconditionVertexLinksInternal() override {
      if(getDimensionality() == 2) {
        preconditionEdgesInternal();
        hasPreconditionedVertexLinks_ = true;
      } else if(getDimensionality() == 3) {
        preconditionTrianglesInternal();
        hasPreconditionedVertexLinks_ = true;
      } else {
        this->printErr("Unsupported dimension for vertex link precondtion");
        return -1;
      }
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() override {
      hasPreconditionedVertexNeighbors_ = true;
      return 0;
    }

    inline int preconditionVertexStarsInternal() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!cellArray_)
        return -1;
#endif

      hasPreconditionedVertexStars_ = true;
      return 0;
    }

    inline int preconditionVertexTrianglesInternal() override {
      preconditionTriangles();
      hasPreconditionedVertexTriangles_ = true;
      return 0;
    }

    /**
     * Initialize the cache with the ratio.
     */
    void initCache(const float ratio = 0.2) {
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
    size_t getCacheSize() {
      return cacheSize_;
    }

    /**
     * Reset the cache size to better fit in the parallel or sequential
     * algorithms.
     */
    int resetCache(int option) {
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
    int clear();

    /**
     * Find the corresponding node index given the id.
     */
    SimplexId findNodeIndex(SimplexId id, int idType) const {
      const std::vector<SimplexId> *intervals = nullptr;
      // determine which vector to search
      if(idType == EDGE_ID) {
        intervals = &edgeIntervals_;
      } else if(idType == TRIANGLE_ID) {
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
    ImplicitCluster *searchCache(const SimplexId &nodeId,
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
                             bool computeInternalEdgeMap) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId edgeCount = 0,
                verticesPerCell = cellArray_->getCellVertexNumber(0);

      if(!nodePtr->internalEdgeMap_.empty()) {
        // if the edge map has been computed before and only request the edge
        // list
        if(computeInternalEdgeList) {
          nodePtr->internalEdgeList_ = std::vector<std::array<SimplexId, 2>>(
            nodePtr->internalEdgeMap_.size());
          for(auto iter = nodePtr->internalEdgeMap_.begin();
              iter != nodePtr->internalEdgeMap_.end(); iter++) {
            nodePtr->internalEdgeList_[iter->second - 1] = iter->first;
          }
          return 0;
        }
      }

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 2> edgeIds;

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          // the edge does not belong to the current node
          if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[1] = cellArray_->getCellVertex(cid, k);

            // not found in the edge map - assign new edge id
            if(nodePtr->internalEdgeMap_.find(edgeIds)
               == nodePtr->internalEdgeMap_.end()) {
              edgeCount++;
              nodePtr->internalEdgeMap_[edgeIds] = edgeCount;
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 2> edgeIds;

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cid, j);
            edgeIds[1] = cellArray_->getCellVertex(cid, k);

            // the edge is in the current node
            if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
               && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
              if(nodePtr->internalEdgeMap_.find(edgeIds)
                 == nodePtr->internalEdgeMap_.end()) {
                edgeCount++;
                (nodePtr->internalEdgeMap_)[edgeIds] = edgeCount;
              }
            }
          }
        }
      }

      if(computeInternalEdgeList) {
        nodePtr->internalEdgeList_ = std::vector<std::array<SimplexId, 2>>(
          nodePtr->internalEdgeMap_.size());
        for(auto iter = nodePtr->internalEdgeMap_.begin();
            iter != nodePtr->internalEdgeMap_.end(); iter++) {
          nodePtr->internalEdgeList_[iter->second - 1] = iter->first;
        }
      }

      if(!computeInternalEdgeMap) {
        nodePtr->internalEdgeMap_.clear();
      }

      return 0;
    }

    /**
     * Build the external edge list in the node.
     */
    int buildExternalEdgeMap(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      boost::unordered_map<SimplexId, std::vector<std::array<SimplexId, 2>>>
        edgeNodes;

      // loop through the external cell list
      for(size_t i = 0; i < externalCells_[nodePtr->nid].size(); i++) {
        std::array<SimplexId, 2> edgeIds;
        SimplexId cellId = externalCells_[nodePtr->nid][i];

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cellId, j);
            edgeIds[1] = cellArray_->getCellVertex(cellId, k);

            // check if the edge is an external edge
            if(edgeIds[0] <= vertexIntervals_[nodePtr->nid - 1]
               && edgeIds[1] > vertexIntervals_[nodePtr->nid - 1]
               && edgeIds[1] <= vertexIntervals_[nodePtr->nid]) {
              SimplexId nid = vertexIndices_[edgeIds[0]];
              edgeNodes[nid].push_back(edgeIds);
            }
          }
        }
      }

      boost::unordered_map<
        SimplexId, std::vector<std::array<SimplexId, 2>>>::iterator iter;
      for(iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(exnode) {
          if(exnode->internalEdgeMap_.empty())
            buildInternalEdgeMap(exnode, false, true);
          for(std::array<SimplexId, 2> edgePair : iter->second) {
            (nodePtr->externalEdgeMap_)[edgePair]
              = exnode->internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          ImplicitCluster tmpCluster(iter->first);
          buildInternalEdgeMap(&tmpCluster, false, true);
          for(std::array<SimplexId, 2> edgePair : iter->second) {
            (nodePtr->externalEdgeMap_)[edgePair]
              = tmpCluster.internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Build the internal triangle list in the node.
     */
    int buildInternalTriangleMap(ImplicitCluster *const nodePtr,
                                 bool computeInternalTriangleList,
                                 bool computeInternalTriangleMap) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId triangleCount = 0, verticesPerCell = 4;
      std::vector<std::vector<SimplexId>> localTriangleStars;

      if(!nodePtr->internalTriangleMap_.empty()) {
        if(computeInternalTriangleList) {
          nodePtr->internalTriangleList_
            = std::vector<std::array<SimplexId, 3>>(
              nodePtr->internalTriangleMap_.size());
          for(auto iter = nodePtr->internalTriangleMap_.begin();
              iter != nodePtr->internalTriangleMap_.end(); iter++) {
            (nodePtr->internalTriangleList_)[iter->second - 1] = iter->first;
          }
          return 0;
        }
      }

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 3> triangleIds;

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          // the triangle does not belong to the current node
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
            for(SimplexId l = k + 1; l < verticesPerCell; l++) {
              triangleIds[1] = cellArray_->getCellVertex(cid, k);
              triangleIds[2] = cellArray_->getCellVertex(cid, l);

              if(nodePtr->internalTriangleMap_.find(triangleIds)
                 == nodePtr->internalTriangleMap_.end()) {
                triangleCount++;
                (nodePtr->internalTriangleMap_)[triangleIds] = triangleCount;
              }
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 3> triangleIds;

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
             && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);

                if(nodePtr->internalTriangleMap_.find(triangleIds)
                   == nodePtr->internalTriangleMap_.end()) {
                  triangleCount++;
                  (nodePtr->internalTriangleMap_)[triangleIds] = triangleCount;
                }
              }
            }
          }
        }
      }

      if(computeInternalTriangleList) {
        nodePtr->internalTriangleList_ = std::vector<std::array<SimplexId, 3>>(
          nodePtr->internalTriangleMap_.size());
        for(auto iter = nodePtr->internalTriangleMap_.begin();
            iter != nodePtr->internalTriangleMap_.end(); iter++) {
          (nodePtr->internalTriangleList_)[iter->second - 1] = iter->first;
        }
      }

      if(!computeInternalTriangleMap) {
        nodePtr->internalTriangleMap_.clear();
      }

      return 0;
    }

    /**
     * Build the external triangle list in the node.
     */
    int buildExternalTriangleMap(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      boost::unordered_map<SimplexId, std::vector<std::array<SimplexId, 3>>>
        nodeTriangles;

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 3> triangleIds;

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] <= vertexIntervals_[nodePtr->nid - 1]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);

                if(triangleIds[1] > vertexIntervals_[nodePtr->nid - 1]
                   && triangleIds[1] <= vertexIntervals_[nodePtr->nid]) {
                  SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                  nodeTriangles[nodeNum].push_back(triangleIds);
                } else if(triangleIds[2] > vertexIntervals_[nodePtr->nid - 1]
                          && triangleIds[2] <= vertexIntervals_[nodePtr->nid]) {
                  SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                  nodeTriangles[nodeNum].push_back(triangleIds);
                }
              }
            }
          }
        }
      }

      boost::unordered_map<
        SimplexId, std::vector<std::array<SimplexId, 3>>>::iterator iter;
      for(iter = nodeTriangles.begin(); iter != nodeTriangles.end(); iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(exnode) {
          if(exnode->internalTriangleMap_.empty())
            buildInternalTriangleMap(exnode, false, true);
          for(std::array<SimplexId, 3> triangleVec : iter->second) {
            (nodePtr->externalTriangleMap_)[triangleVec]
              = exnode->internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[iter->first - 1];
          }
        } else {
          ImplicitCluster tmpCluster(iter->first);
          buildInternalTriangleMap(&tmpCluster, false, true);
          for(std::array<SimplexId, 3> triangleVec : iter->second) {
            (nodePtr->externalTriangleMap_)[triangleVec]
              = tmpCluster.internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the number of internal edges in the node.
     */
    SimplexId countInternalEdges(SimplexId nodeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodeId <= 0 || nodeId > nodeNumber_)
        return -1;
#endif

      SimplexId edgeCount = 0,
                verticesPerCell = cellArray_->getCellVertexNumber(0);
      boost::unordered_set<std::array<SimplexId, 2>,
                           boost::hash<std::array<SimplexId, 2>>>
        edgeSet;

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        std::array<SimplexId, 2> edgeIds;

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          // the edge does not belong to the current node
          if(edgeIds[0] > vertexIntervals_[nodeId]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[1] = cellArray_->getCellVertex(cid, k);

            // not found in the edge map - assign new edge id
            if(edgeSet.find(edgeIds) == edgeSet.end()) {
              edgeCount++;
              edgeSet.insert(edgeIds);
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodeId]) {
        std::array<SimplexId, 2> edgeIds;

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cid, j);
            edgeIds[1] = cellArray_->getCellVertex(cid, k);

            // the edge is in the current node
            if(edgeIds[0] > vertexIntervals_[nodeId - 1]
               && edgeIds[0] <= vertexIntervals_[nodeId]) {
              if(edgeSet.find(edgeIds) == edgeSet.end()) {
                edgeCount++;
                edgeSet.insert(edgeIds);
              }
            }
          }
        }
      }

      return edgeCount;
    }

    /**
     * Get the number of internal triangles in the node.
     */
    int countInternalTriangles(SimplexId nodeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodeId <= 0 || nodeId > nodeNumber_)
        return -1;
#endif

      SimplexId triangleCount = 0,
                verticesPerCell = cellArray_->getCellVertexNumber(0);
      boost::unordered_set<std::vector<SimplexId>,
                           boost::hash<std::vector<SimplexId>>>
        triangleSet;

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        std::vector<SimplexId> triangleIds(3);

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          // the triangle does not belong to the current node
          if(triangleIds[0] > vertexIntervals_[nodeId]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
            for(SimplexId l = k + 1; l < verticesPerCell; l++) {
              triangleIds[1] = cellArray_->getCellVertex(cid, k);
              triangleIds[2] = cellArray_->getCellVertex(cid, l);

              if(triangleSet.find(triangleIds) == triangleSet.end()) {
                triangleCount++;
                triangleSet.insert(triangleIds);
              }
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodeId]) {
        std::vector<SimplexId> triangleIds(3);

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodeId - 1]
             && triangleIds[0] <= vertexIntervals_[nodeId]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);

                if(triangleSet.find(triangleIds) == triangleSet.end()) {
                  triangleCount++;
                  triangleSet.insert(triangleIds);
                }
              }
            }
          }
        }
      }

      return triangleCount;
    }

    /**
     * Get the cell edges for all cells in a given cluster.
     * Check if the tetraEdges_ is NULL before calling the function.
     */
    int getClusterCellEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = 4;
      nodePtr->tetraEdges_ = std::vector<std::array<SimplexId, 6>>(
        cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);
      boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
        edgeNodes;

      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }

      for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
          i <= cellIntervals_[nodePtr->nid]; i++) {
        int cnt = 0;
        // get the internal edge id from the map
        for(SimplexId k = 1; k < verticesPerCell; k++) {
          std::array<SimplexId, 2> edgePair
            = {(SimplexId)cellArray_->getCellVertex(i, 0),
               (SimplexId)cellArray_->getCellVertex(i, k)};
          (nodePtr
             ->tetraEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][cnt++]
            = nodePtr->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
        }
        for(SimplexId j = 1; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            std::array<SimplexId, 2> edgePair
              = {(SimplexId)cellArray_->getCellVertex(i, j),
                 (SimplexId)cellArray_->getCellVertex(i, k)};
            if(edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
              (nodePtr
                 ->tetraEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][cnt++]
                = nodePtr->internalEdgeMap_.at(edgePair)
                  + edgeIntervals_[nodePtr->nid - 1];
            }
            // group the external edges by node id
            else {
              std::vector<SimplexId> edgeTuple{
                i, cnt++, edgePair[0], edgePair[1]};
              SimplexId nodeNum = vertexIndices_[edgePair[0]];
              edgeNodes[nodeNum].push_back(edgeTuple);
            }
          }
        }
      }

      for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(exnode) {
          if(exnode->internalEdgeMap_.empty())
            buildInternalEdgeMap(exnode, false, true);
          for(std::vector<SimplexId> edgeTuple : iter->second) {
            std::array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
            (nodePtr
               ->tetraEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                              - 1][edgeTuple[1]]
              = exnode->internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          ImplicitCluster tmpCluster(iter->first);
          buildInternalEdgeMap(&tmpCluster, false, true);
          for(std::vector<SimplexId> edgeTuple : iter->second) {
            std::array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
            (nodePtr
               ->tetraEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                              - 1][edgeTuple[1]]
              = tmpCluster.internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the cell neighbors for all cells in a given cluster.
     */
    int getClusterCellNeighbors(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      std::vector<std::vector<SimplexId>> localCellNeighbors(
        cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);
      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      std::vector<std::vector<SimplexId>> localVertexStars;

      if(nodePtr->vertexStars_.empty()) {
        getClusterVertexStars(nodePtr);
      }
      // sort the vertex star vector
      nodePtr->vertexStars_.copyTo(localVertexStars);
      for(size_t i = 0; i < localVertexStars.size(); i++) {
        sort(localVertexStars[i].begin(), localVertexStars[i].end());
      }

      boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
        nodeMaps;
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        for(SimplexId j = 1; j < verticesPerCell; j++) {
          if(cellArray_->getCellVertex(cid, j)
             > vertexIntervals_[nodePtr->nid]) {
            SimplexId nodeId
              = vertexIndices_[cellArray_->getCellVertex(cid, j)];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              ImplicitCluster newNode(nodeId);
              getClusterVertexStars(&newNode);
              // sort the vertex star list
              std::vector<std::vector<SimplexId>> tmpVec;
              newNode.vertexStars_.copyTo(tmpVec);
              for(size_t i = 0; i < tmpVec.size(); i++) {
                sort(tmpVec[i].begin(), tmpVec[i].end());
              }
              nodeMaps[nodeId] = tmpVec;
            }
          }
        }
      }

      if(getDimensionality() == 2) {
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          for(SimplexId j = 0; j < 3; j++) {

            SimplexId v0 = cellArray_->getCellVertex(cid, j);
            SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 3);
            SimplexId localV0 = v0 - vertexIntervals_[nodePtr->nid - 1] - 1;
            SimplexId localV1 = v1 - vertexIntervals_[nodePtr->nid - 1] - 1;

            std::vector<SimplexId> stars0, stars1;
            if(v0 <= vertexIntervals_[nodePtr->nid]) {
              stars0 = localVertexStars[localV0];
            } else {
              localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
              stars0 = nodeMaps[vertexIndices_[v0]][localV0];
            }
            if(v1 <= vertexIntervals_[nodePtr->nid]) {
              stars1 = localVertexStars[localV1];
            } else {
              localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
              stars1 = nodeMaps[vertexIndices_[v1]][localV1];
            }

            // perform an intersection of the 2 sorted star lists
            SimplexId pos0 = 0, pos1 = 0;
            SimplexId intersection = -1;

            while((pos0 < (SimplexId)stars0.size())
                  && (pos1 < (SimplexId)stars1.size())) {
              SimplexId biggest = stars0[pos0];
              if(stars1[pos1] > biggest) {
                biggest = stars1[pos1];
              }

              for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
                if(stars0[l] < biggest) {
                  pos0++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
                if(stars1[l] < biggest) {
                  pos1++;
                } else {
                  break;
                }
              }

              if(pos0 >= (SimplexId)stars0.size()
                 || pos1 >= (SimplexId)stars1.size())
                break;

              if(stars0[pos0] == stars1[pos1]) {
                if(stars0[pos0] != cid) {
                  intersection = stars0[pos0];
                  break;
                }
                pos0++;
                pos1++;
              }
            }

            if(intersection != -1) {
              localCellNeighbors[cid - cellIntervals_[nodePtr->nid - 1] - 1]
                .push_back(intersection);
            }
          }
        }
      }

      else if(getDimensionality() == 3) {
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          // go triangle by triangle
          for(SimplexId j = 0; j < 4; j++) {

            SimplexId v0 = cellArray_->getCellVertex(cid, j % 4);
            SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 4);
            SimplexId v2 = cellArray_->getCellVertex(cid, (j + 2) % 4);

            SimplexId localV0 = v0 - vertexIntervals_[nodePtr->nid - 1] - 1;
            SimplexId localV1 = v1 - vertexIntervals_[nodePtr->nid - 1] - 1;
            SimplexId localV2 = v2 - vertexIntervals_[nodePtr->nid - 1] - 1;

            std::vector<SimplexId> stars0, stars1, stars2;
            if(v0 <= vertexIntervals_[nodePtr->nid]) {
              stars0 = localVertexStars[localV0];
            } else {
              localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
              stars0 = nodeMaps[vertexIndices_[v0]][localV0];
            }
            if(v1 <= vertexIntervals_[nodePtr->nid]) {
              stars1 = localVertexStars[localV1];
            } else {
              localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
              stars1 = nodeMaps[vertexIndices_[v1]][localV1];
            }
            if(v2 <= vertexIntervals_[nodePtr->nid]) {
              stars2 = localVertexStars[localV2];
            } else {
              localV2 = v2 - vertexIntervals_[vertexIndices_[v2] - 1] - 1;
              stars2 = nodeMaps[vertexIndices_[v2]][localV2];
            }

            // perform an intersection of the 3 (sorted) star lists
            SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
            SimplexId intersection = -1;

            while((pos0 < (SimplexId)stars0.size())
                  && (pos1 < (SimplexId)stars1.size())
                  && (pos2 < (SimplexId)stars2.size())) {

              SimplexId biggest = stars0[pos0];
              if(stars1[pos1] > biggest) {
                biggest = stars1[pos1];
              }
              if(stars2[pos2] > biggest) {
                biggest = stars2[pos2];
              }

              for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
                if(stars0[l] < biggest) {
                  pos0++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
                if(stars1[l] < biggest) {
                  pos1++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos2; l < (SimplexId)stars2.size(); l++) {
                if(stars2[l] < biggest) {
                  pos2++;
                } else {
                  break;
                }
              }

              if(pos0 >= (SimplexId)stars0.size()
                 || pos1 >= (SimplexId)stars1.size()
                 || pos2 >= (SimplexId)stars2.size())
                break;

              if((stars0[pos0] == stars1[pos1])
                 && (stars0[pos0] == stars2[pos2])) {
                if(stars0[pos0] != cid) {
                  intersection = stars0[pos0];
                  break;
                }
                pos0++;
                pos1++;
                pos2++;
              }
            }

            if(intersection != -1) {
              localCellNeighbors[cid - cellIntervals_[nodePtr->nid - 1] - 1]
                .push_back(intersection);
            }
          }
        }
      }

      nodePtr->cellNeighbors_.fillFrom(localCellNeighbors);
      return 0;
    }

    /**
     * Get the cell triangles for all cells in a given cluster.
     */
    int getClusterCellTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = 4;
      boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
        nodeTriangles;

      nodePtr->tetraTriangles_ = std::vector<std::array<SimplexId, 4>>(
        cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);

      if(nodePtr->internalTriangleMap_.empty()) {
        buildInternalTriangleMap(nodePtr, false, true);
      }

      for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
          i <= cellIntervals_[nodePtr->nid]; i++) {
        std::array<SimplexId, 3> triangleVec;
        // get the internal triangle from the map
        triangleVec[0] = cellArray_->getCellVertex(i, 0);
        for(SimplexId k = 1; k < verticesPerCell - 1; k++) {
          triangleVec[1] = cellArray_->getCellVertex(i, k);
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleVec[2] = cellArray_->getCellVertex(i, l);
            (nodePtr->tetraTriangles_)[i - cellIntervals_[nodePtr->nid - 1] - 1]
                                      [k + l - 3]
              = nodePtr->internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[nodePtr->nid - 1];
          }
        }
        // group the external triangles by node id
        triangleVec[0] = cellArray_->getCellVertex(i, 1);
        triangleVec[1] = cellArray_->getCellVertex(i, 2);
        triangleVec[2] = cellArray_->getCellVertex(i, 3);
        if(triangleVec[0] <= vertexIntervals_[nodePtr->nid]) {
          (nodePtr->tetraTriangles_)[i - cellIntervals_[nodePtr->nid - 1] - 1]
            .back()
            = nodePtr->internalTriangleMap_.at(triangleVec)
              + triangleIntervals_[nodePtr->nid - 1];
        } else {
          std::vector<SimplexId> triangleTuple
            = {i, triangleVec[0], triangleVec[1], triangleVec[2]};
          SimplexId nodeNum = vertexIndices_[triangleVec[0]];
          nodeTriangles[nodeNum].push_back(triangleTuple);
        }
      }

      for(auto iter = nodeTriangles.begin(); iter != nodeTriangles.end();
          iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(exnode) {
          if(exnode->internalTriangleMap_.empty())
            buildInternalTriangleMap(exnode, false, true);
          for(std::vector<SimplexId> triangleVec : iter->second) {
            std::array<SimplexId, 3> triangle
              = {triangleVec[1], triangleVec[2], triangleVec[3]};
            (nodePtr->tetraTriangles_)[triangleVec[0]
                                       - cellIntervals_[nodePtr->nid - 1] - 1]
              .back()
              = exnode->internalTriangleMap_.at(triangle)
                + triangleIntervals_[iter->first - 1];
          }
        } else {
          ImplicitCluster tmpCluster(iter->first);
          buildInternalTriangleMap(&tmpCluster, false, true);
          for(std::vector<SimplexId> triangleVec : iter->second) {
            std::array<SimplexId, 3> triangle
              = {triangleVec[1], triangleVec[2], triangleVec[3]};
            (nodePtr->tetraTriangles_)[triangleVec[0]
                                       - cellIntervals_[nodePtr->nid - 1] - 1]
              .back()
              = tmpCluster.internalTriangleMap_.at(triangle)
                + triangleIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the edge links for all the edges in a given cluster.
     */
    int getClusterEdgeLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localEdgeNum
        = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
        linksCount(localEdgeNum, 0);

      if(getDimensionality() == 2) {
        if(nodePtr->edgeStars_.empty()) {
          getClusterEdgeStars(nodePtr);
        }
        // set the offsets vector
        boost::unordered_map<std::array<SimplexId, 2>,
                             SimplexId>::const_iterator iter;
        for(iter = nodePtr->internalEdgeMap_.begin();
            iter != nodePtr->internalEdgeMap_.end(); iter++) {
          for(SimplexId j = 0; j < nodePtr->edgeStars_.size(iter->second - 1);
              j++) {
            SimplexId cellId = nodePtr->edgeStars_.get(iter->second - 1, j);
            for(int k = 0; k < 3; k++) {
              SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
              if((vertexId != iter->first[0]) && (vertexId != iter->first[1])) {
                offsets[iter->second]++;
                break;
              }
            }
          }
        }

        // compute partial sum of number of links per edge
        for(SimplexId i = 1; i <= localEdgeNum; i++) {
          offsets[i] += offsets[i - 1];
        }

        // allocate the flat vector for edge link data
        std::vector<SimplexId> edgeLinkData(offsets.back());

        // fill the flat vector using offsets and count vectors
        for(iter = nodePtr->internalEdgeMap_.begin();
            iter != nodePtr->internalEdgeMap_.end(); iter++) {
          for(SimplexId j = 0; j < nodePtr->edgeStars_.size(iter->second - 1);
              j++) {
            SimplexId cellId = nodePtr->edgeStars_.get(iter->second - 1, j);
            for(int k = 0; k < 3; k++) {
              SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
              if((vertexId != iter->first[0]) && (vertexId != iter->first[1])) {
                SimplexId localVertexId
                  = vertexId - vertexIntervals_[nodePtr->nid - 1] - 1;
                edgeLinkData[offsets[localVertexId] + linksCount[localVertexId]]
                  = vertexId;
                break;
              }
            }
          }
        }

        // fill FlatJaggedArray struct
        nodePtr->edgeLinks_.setData(
          std::move(edgeLinkData), std::move(offsets));
      } else if(getDimensionality() == 3) {
        if(nodePtr->tetraEdges_.empty()) {
          getClusterCellEdges(nodePtr);
        }
        if(nodePtr->internalEdgeMap_.empty()) {
          buildInternalEdgeMap(nodePtr, false, true);
        }

        // set the offsets vector
        SimplexId localCellNum
          = cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1];
        for(SimplexId cid = 0; cid < localCellNum; cid++) {
          SimplexId cellId = cid + cellIntervals_[nodePtr->nid - 1] + 1;
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
               (SimplexId)cellArray_->getCellVertex(cellId, 1),
               (SimplexId)cellArray_->getCellVertex(cellId, 2),
               (SimplexId)cellArray_->getCellVertex(cellId, 3)};
          std::array<SimplexId, 2> edgePair;
          edgePair[0] = vertexIds[0];
          for(SimplexId j = 1; j < 4; j++) {
            edgePair[1] = vertexIds[j];
            offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
          }
          if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = vertexIds[1];
            for(int j = 2; j < 4; j++) {
              edgePair[1] = vertexIds[j];
              offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
            }
            if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
              edgePair = {vertexIds[2], vertexIds[3]};
              offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          std::array<SimplexId, 2> edgeIds;
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2),
               (SimplexId)cellArray_->getCellVertex(cid, 3)};

          // loop through each edge of the cell
          for(SimplexId j = 0; j < 3; j++) {
            for(SimplexId k = j + 1; k < 4; k++) {
              edgeIds[0] = vertexIds[j];
              edgeIds[1] = vertexIds[k];

              // the edge is in the current node
              if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
                 && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
                offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
              }
            }
          }
        }

        // compute partial sum of number of neighbors per vertex
        for(SimplexId i = 1; i <= localEdgeNum; i++) {
          offsets[i] += offsets[i - 1];
        }

        // allocate the flat vector for edge link data
        std::vector<SimplexId> edgeLinkData(offsets.back());

        // fill the flat vector using offsets and count vectors
        for(SimplexId cid = 0; cid < localCellNum; cid++) {
          SimplexId cellId = cid + cellIntervals_[nodePtr->nid - 1] + 1;
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
               (SimplexId)cellArray_->getCellVertex(cellId, 1),
               (SimplexId)cellArray_->getCellVertex(cellId, 2),
               (SimplexId)cellArray_->getCellVertex(cellId, 3)};
          std::array<SimplexId, 2> edgePair;
          edgePair[0] = vertexIds[0];
          for(SimplexId j = 1; j < 4; j++) {
            edgePair[1] = vertexIds[j];
            SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgePair) - 1;
            edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
              = nodePtr->tetraEdges_.at(cid)[6 - j];
            linksCount[localEdgeId]++;
          }
          if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = vertexIds[1];
            for(int j = 2; j < 4; j++) {
              edgePair[1] = vertexIds[j];
              SimplexId localEdgeId
                = nodePtr->internalEdgeMap_.at(edgePair) - 1;
              edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
                = nodePtr->tetraEdges_.at(cid)[4 - j];
              linksCount[localEdgeId]++;
            }
            if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
              edgePair = {vertexIds[2], vertexIds[3]};
              SimplexId localEdgeId
                = nodePtr->internalEdgeMap_.at(edgePair) - 1;
              edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
                = nodePtr->tetraEdges_.at(cid)[0];
              linksCount[localEdgeId]++;
            }
          }
        }

        // loop through the external cell list
        boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2),
               (SimplexId)cellArray_->getCellVertex(cid, 3)};

          std::array<SimplexId, 2> edgeIds;
          // loop through each edge of the cell
          for(SimplexId j = 0; j < 3; j++) {
            for(SimplexId k = j + 1; k < 4; k++) {
              edgeIds[0] = vertexIds[j];
              edgeIds[1] = vertexIds[k];

              // the edge is in the current node
              if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
                 && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
                std::array<SimplexId, 2> otherEdge = {-1, -1};
                for(int i = 0; i < 4; i++) {
                  if(vertexIds[i] != edgeIds[0] && vertexIds[i] != edgeIds[1]) {
                    if(otherEdge[0] == -1) {
                      otherEdge[0] = vertexIds[i];
                    } else if(otherEdge[1] == -1) {
                      otherEdge[1] = vertexIds[i];
                    } else {
                      printErr("[CompactTriangulation] More than two other "
                               "vertices are "
                               "found in the edge!\n");
                    }
                  }
                }
                SimplexId nodeId = vertexIndices_[otherEdge[0]];
                if(nodeMaps.find(nodeId) == nodeMaps.end()) {
                  nodeMaps[nodeId] = ImplicitCluster(nodeId);
                  buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
                }
                SimplexId localEdgeId
                  = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
                edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
                  = nodeMaps[nodeId].internalEdgeMap_.at(otherEdge);
                linksCount[localEdgeId]++;
              }
            }
          }
        }

        // fill FlatJaggedArray struct
        nodePtr->edgeLinks_.setData(
          std::move(edgeLinkData), std::move(offsets));
      }

      return 0;
    }

    /**
     * Get the edge stars for all the edges in a given cluster.
     */
    int getClusterEdgeStars(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localEdgeNum
        = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
        starsCount(localEdgeNum);

      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }

      // set the offsets vector
      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 2> edgeIds;
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          // the edge does not belong to the current node
          if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[1] = cellArray_->getCellVertex(cid, k);
            offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
          }
        }
      }
      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 2> edgeIds;
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cid, j);
            edgeIds[1] = cellArray_->getCellVertex(cid, k);
            if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
               && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
              offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
            }
          }
        }
      }

      // compute partial sum of number of stars per edge
      for(SimplexId i = 1; i <= localEdgeNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for edge star data
      std::vector<SimplexId> edgeStarData(offsets.back());

      // fill the flat vector using offsets and count vectors
      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 2> edgeIds;
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          // the edge does not belong to the current node
          if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[1] = cellArray_->getCellVertex(cid, k);
            SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
            edgeStarData[offsets[localEdgeId] + starsCount[localEdgeId]] = cid;
            starsCount[localEdgeId]++;
          }
        }
      }
      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 2> edgeIds;
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cid, j);
            edgeIds[1] = cellArray_->getCellVertex(cid, k);
            if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
               && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
              SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
              edgeStarData[offsets[localEdgeId] + starsCount[localEdgeId]]
                = cid;
              starsCount[localEdgeId]++;
            }
          }
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->edgeStars_.setData(std::move(edgeStarData), std::move(offsets));

      return 0;
    }

    /**
     * Get the edge triangles for all the edges in a given cluster.
     */
    int getClusterEdgeTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localEdgeNum
        = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
        trianglesCount(localEdgeNum, 0);

      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }
      if(nodePtr->internalTriangleMap_.empty()) {
        buildInternalTriangleMap(nodePtr, false, true);
      }
      if(nodePtr->externalTriangleMap_.empty()) {
        buildExternalTriangleMap(nodePtr);
      }

      // set the offsets vector
      boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        std::array<SimplexId, 2> edge1 = {iter->first[0], iter->first[1]};
        std::array<SimplexId, 2> edge2 = {iter->first[0], iter->first[2]};
        offsets[nodePtr->internalEdgeMap_.at(edge1)]++;
        offsets[nodePtr->internalEdgeMap_.at(edge2)]++;
        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          edge1 = {iter->first[1], iter->first[2]};
          offsets[nodePtr->internalEdgeMap_.at(edge1)]++;
        }
      }
      for(iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        std::array<SimplexId, 2> edge = {iter->first.at(1), iter->first.at(2)};
        if(edge[0] > vertexIntervals_[nodePtr->nid - 1]
           && edge[0] <= vertexIntervals_[nodePtr->nid]) {
          offsets[nodePtr->internalEdgeMap_.at(edge)]++;
        }
      }

      // compute partial sum of number of triangles per edge
      for(SimplexId i = 1; i <= localEdgeNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for edge triangle data
      std::vector<SimplexId> edgeTriangleData(offsets.back());

      // fill the flat vector using offsets and count vectors
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        std::array<SimplexId, 2> edge1 = {iter->first[0], iter->first[1]};
        std::array<SimplexId, 2> edge2 = {iter->first[0], iter->first[2]};
        SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edge1) - 1;
        edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
          = iter->second + triangleIntervals_[nodePtr->nid - 1];
        trianglesCount[localEdgeId]++;
        localEdgeId = nodePtr->internalEdgeMap_.at(edge2) - 1;
        edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
          = iter->second + triangleIntervals_[nodePtr->nid - 1];
        trianglesCount[localEdgeId]++;

        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          edge1 = {iter->first[1], iter->first[2]};
          localEdgeId = nodePtr->internalEdgeMap_.at(edge1) - 1;
          edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
            = iter->second + triangleIntervals_[nodePtr->nid - 1];
          trianglesCount[localEdgeId]++;
        }
      }

      // for external triangles
      for(iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        std::array<SimplexId, 2> edge = {iter->first.at(1), iter->first.at(2)};
        if(edge[0] > vertexIntervals_[nodePtr->nid - 1]
           && edge[0] <= vertexIntervals_[nodePtr->nid]) {
          SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edge) - 1;
          edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
            = iter->second;
          trianglesCount[localEdgeId]++;
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->edgeTriangles_.setData(
        std::move(edgeTriangleData), std::move(offsets));

      return 0;
    }

    /**
     * Get the triangle edges for all the triangles in a given cluster.
     */
    int getClusterTriangleEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->triangleEdges_ = std::vector<std::array<SimplexId, 3>>(
        triangleIntervals_[nodePtr->nid]
        - triangleIntervals_[nodePtr->nid - 1]);
      boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
        edgeNodes;

      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }
      if(nodePtr->internalTriangleMap_.empty()) {
        buildInternalTriangleMap(nodePtr, false, true);
      }

      for(auto iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        // since the first vertex of the triangle is in the node ...
        std::array<SimplexId, 2> edgePair = {iter->first[0], iter->first[1]};
        (nodePtr->triangleEdges_)[iter->second - 1][0]
          = nodePtr->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
        edgePair[1] = iter->first[2];
        (nodePtr->triangleEdges_)[iter->second - 1][1]
          = nodePtr->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
        edgePair[0] = iter->first[1];
        if(edgePair[0] > vertexIntervals_[nodePtr->nid - 1]
           && edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
          (nodePtr->triangleEdges_)[iter->second - 1][2]
            = nodePtr->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
        } else {
          std::vector<SimplexId> edgeTuple{
            iter->second - 1, edgePair[0], edgePair[1]};
          SimplexId nodeNum = vertexIndices_[edgePair[0]];
          edgeNodes[nodeNum].push_back(edgeTuple);
        }
      }

      for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(exnode) {
          if(exnode->internalEdgeMap_.empty())
            buildInternalEdgeMap(exnode, false, true);
          for(std::vector<SimplexId> edgeTuple : iter->second) {
            std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
            (nodePtr->triangleEdges_)[edgeTuple[0]][2]
              = exnode->internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          ImplicitCluster tmpCluster(iter->first);
          buildInternalEdgeMap(&tmpCluster, false, true);
          for(std::vector<SimplexId> edgeTuple : iter->second) {
            std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
            (nodePtr->triangleEdges_)[edgeTuple[0]][2]
              = tmpCluster.internalEdgeMap_.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the triangle links for all the triangles in a given cluster.
     */
    int getClusterTriangleLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localTriangleNum = triangleIntervals_[nodePtr->nid]
                                   - triangleIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localTriangleNum + 1, 0),
        linksCount(localTriangleNum, 0);

      if(nodePtr->triangleStars_.empty()) {
        getClusterTriangleStars(nodePtr);
      }

      // set the offsets vector
      boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::const_iterator
        iter;
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        for(SimplexId i = 0; i < nodePtr->triangleStars_.size(iter->second - 1);
            i++) {
          SimplexId cellId = nodePtr->triangleStars_.get(iter->second - 1, i);
          for(int j = 0; j < 4; j++) {
            SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
            if((vertexId != iter->first[0]) && (vertexId != iter->first[1])
               && (vertexId != iter->first[2])) {
              offsets[iter->second]++;
              break;
            }
          }
        }
      }

      // compute partial sum of number of links per triangle
      for(SimplexId i = 1; i <= localTriangleNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for triangle link data
      std::vector<SimplexId> triangleLinkData(offsets.back());

      // fill the flat vector using offsets and count vectors
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        for(SimplexId i = 0; i < nodePtr->triangleStars_.size(iter->second - 1);
            i++) {
          SimplexId cellId = nodePtr->triangleStars_.get(iter->second - 1, i);
          for(int j = 0; j < 4; j++) {
            SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
            if((vertexId != iter->first[0]) && (vertexId != iter->first[1])
               && (vertexId != iter->first[2])) {
              triangleLinkData[offsets[iter->second - 1]
                               + linksCount[iter->second - 1]]
                = vertexId;
              linksCount[iter->second - 1]++;
              break;
            }
          }
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->triangleLinks_.setData(
        std::move(triangleLinkData), std::move(offsets));

      return 0;
    }

    /**
     * Get the triangle stars for all the triangles in a given cluster.
     */
    int getClusterTriangleStars(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localTriangleNum = triangleIntervals_[nodePtr->nid]
                                   - triangleIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localTriangleNum + 1, 0),
        starsCount(localTriangleNum, 0);

      if(nodePtr->internalTriangleMap_.empty()) {
        buildInternalTriangleMap(nodePtr, false, true);
      }

      // set the offsets vector
      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 3> triangleIds;
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
            for(SimplexId l = k + 1; l < verticesPerCell; l++) {
              triangleIds[1] = cellArray_->getCellVertex(cid, k);
              triangleIds[2] = cellArray_->getCellVertex(cid, l);
              offsets[nodePtr->internalTriangleMap_.at(triangleIds)]++;
            }
          }
        }
      }
      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 3> triangleIds;
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
             && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);
                offsets[nodePtr->internalTriangleMap_.at(triangleIds)]++;
              }
            }
          }
        }
      }

      // compute partial sum of number of stars per triangle
      for(SimplexId i = 1; i <= localTriangleNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for triangle star data
      std::vector<SimplexId> triangleStarData(offsets.back());

      // fill the flat vector using offsets and count vectors
      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        std::array<SimplexId, 3> triangleIds;
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          // the triangle does not belong to the current node
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
            break;
          }
          for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
            for(SimplexId l = k + 1; l < verticesPerCell; l++) {
              triangleIds[1] = cellArray_->getCellVertex(cid, k);
              triangleIds[2] = cellArray_->getCellVertex(cid, l);
              SimplexId localTriangleId
                = nodePtr->internalTriangleMap_.at(triangleIds) - 1;
              triangleStarData[offsets[localTriangleId]
                               + starsCount[localTriangleId]]
                = cid;
              starsCount[localTriangleId]++;
            }
          }
        }
      }
      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        std::array<SimplexId, 3> triangleIds;

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
             && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);
                SimplexId localTriangleId
                  = nodePtr->internalTriangleMap_.at(triangleIds) - 1;
                triangleStarData[offsets[localTriangleId]
                                 + starsCount[localTriangleId]]
                  = cid;
                starsCount[localTriangleId]++;
              }
            }
          }
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->triangleStars_.setData(
        std::move(triangleStarData), std::move(offsets));

      return 0;
    }

    /**
     * Get the vertex edges for all the vertices in a given cluster.
     */
    int getClusterVertexEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localVertexNum + 1, 0),
        edgesCount(localVertexNum, 0);

      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }
      if(nodePtr->externalEdgeMap_.empty()) {
        buildExternalEdgeMap(nodePtr);
      }

      // set the offsets vector
      boost::unordered_map<std::array<SimplexId, 2>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalEdgeMap_.begin();
          iter != nodePtr->internalEdgeMap_.end(); iter++) {
        offsets[iter->first[0] - vertexIntervals_[nodePtr->nid - 1]]++;
        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          offsets[iter->first[1] - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
      for(iter = nodePtr->externalEdgeMap_.begin();
          iter != nodePtr->externalEdgeMap_.end(); iter++) {
        offsets[iter->first[1] - vertexIntervals_[nodePtr->nid - 1]]++;
      }
      // compute partial sum of number of edges per vertex
      for(SimplexId i = 1; i <= localVertexNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for vertex edge data
      std::vector<SimplexId> vertexEdgeData(offsets.back());

      // fill the flat vector using offsets and count vectors
      for(iter = nodePtr->internalEdgeMap_.begin();
          iter != nodePtr->internalEdgeMap_.end(); iter++) {
        SimplexId localVertexId
          = iter->first[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
          = iter->second + edgeIntervals_[nodePtr->nid - 1];
        edgesCount[localVertexId]++;
        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          localVertexId
            = iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
          vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
            = iter->second + edgeIntervals_[nodePtr->nid - 1];
          edgesCount[localVertexId]++;
        }
      }
      for(iter = nodePtr->externalEdgeMap_.begin();
          iter != nodePtr->externalEdgeMap_.end(); iter++) {
        SimplexId localVertexId
          = iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
          = iter->second;
        edgesCount[localVertexId]++;
      }

      // fill FlatJaggedArray struct
      nodePtr->vertexEdges_.setData(
        std::move(vertexEdgeData), std::move(offsets));

      return 0;
    }

    /**
     * Get the vertex links for all the vertices in a given cluster.
     */
    int getClusterVertexLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localVertexNum + 1, 0),
        linksCount(localVertexNum, 0);
      boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
      // triangle mesh
      if(getDimensionality() == 2) {
        if(nodePtr->internalEdgeMap_.empty()) {
          buildInternalEdgeMap(nodePtr, false, true);
        }

        // set the offsets vector
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          offsets[cellArray_->getCellVertex(cid, 0)
                  - vertexIntervals_[nodePtr->nid - 1]]++;
          for(SimplexId j = 1; j < 3; j++) {
            if(cellArray_->getCellVertex(cid, j)
               <= vertexIntervals_[nodePtr->nid]) {
              offsets[cellArray_->getCellVertex(cid, j)
                      - vertexIntervals_[nodePtr->nid - 1]]++;
            }
          }
        }
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          for(SimplexId j = 1; j < 3; j++) {
            SimplexId vertexId = cellArray_->getCellVertex(cid, j);
            if(vertexId > vertexIntervals_[nodePtr->nid - 1]
               && vertexId <= vertexIntervals_[nodePtr->nid]) {
              offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
            }
          }
        }

        // compute partial sum of number of links per vertex
        for(SimplexId i = 1; i <= localVertexNum; i++) {
          offsets[i] += offsets[i - 1];
        }

        // allocate the flat vector for vertex link data
        std::vector<SimplexId> vertexLinkData(offsets.back());

        // fill the flat vector using offsets and count vectors
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          std::array<SimplexId, 3> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2)};
          // the first vertex of the cell must be in the cluster
          std::array<SimplexId, 2> edgePair = {vertexIds[1], vertexIds[2]};
          SimplexId nodeId = vertexIndices_[vertexIds[1]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            nodeMaps[nodeId] = ImplicitCluster(nodeId);
            buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
          }
          SimplexId localVertexId
            = vertexIds[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
              + edgeIntervals_[nodeId - 1];
          linksCount[localVertexId]++;

          if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = vertexIds[0];
            localVertexId
              = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodePtr->internalEdgeMap_.at(edgePair)
                + edgeIntervals_[nodePtr->nid - 1];
            linksCount[localVertexId]++;
            if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
              localVertexId
                = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
              edgePair[1] = vertexIds[1];
              vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
                = nodePtr->internalEdgeMap_.at(edgePair)
                  + edgeIntervals_[nodePtr->nid - 1];
              linksCount[localVertexId]++;
            }
          }
        }
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          std::array<SimplexId, 3> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2)};
          std::array<SimplexId, 2> edgePair = {vertexIds[0], vertexIds[2]};
          SimplexId nodeId = vertexIndices_[edgePair[0]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            nodeMaps[nodeId] = ImplicitCluster(nodeId);
            buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
          }
          SimplexId localVertexId
            = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
          if(vertexIds[1] > vertexIntervals_[nodePtr->nid - 1]
             && vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
                + edgeIntervals_[nodeId - 1];
            linksCount[localVertexId]++;
          }
          if(vertexIds[2] > vertexIntervals_[nodePtr->nid - 1]
             && vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
            edgePair[1] = vertexIds[1];
            localVertexId
              = vertexIds[2] - -vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
                + edgeIntervals_[nodeId - 1];
            linksCount[localVertexId]++;
          }
        }

        // fill FlatJaggedArray struct
        nodePtr->vertexLinks_.setData(
          std::move(vertexLinkData), std::move(offsets));

        // tetrahedral mesh
      } else if(getDimensionality() == 3) {

        if(nodePtr->internalTriangleMap_.empty()) {
          buildInternalTriangleMap(nodePtr, false, true);
        }

        // set the offsets vector
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          offsets[cellArray_->getCellVertex(cid, 0)
                  - vertexIntervals_[nodePtr->nid - 1]]++;
          for(SimplexId j = 1; j < 4; j++) {
            if(cellArray_->getCellVertex(cid, j)
               <= vertexIntervals_[nodePtr->nid]) {
              offsets[cellArray_->getCellVertex(cid, j)
                      - vertexIntervals_[nodePtr->nid - 1]]++;
            }
          }
        }
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          for(SimplexId j = 1; j < 4; j++) {
            SimplexId vertexId = cellArray_->getCellVertex(cid, j);
            if(vertexId > vertexIntervals_[nodePtr->nid - 1]
               && vertexId <= vertexIntervals_[nodePtr->nid]) {
              offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
            }
          }
        }

        // compute partial sum of number of links per vertex
        for(SimplexId i = 1; i <= localVertexNum; i++) {
          offsets[i] += offsets[i - 1];
        }

        // allocate the flat vector for vertex link data
        std::vector<SimplexId> vertexLinkData(offsets.back());

        // fill the flat vector using offsets and count vectors
        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2),
               (SimplexId)cellArray_->getCellVertex(cid, 3)};

          // v1: (v2, v3, v4)
          std::array<SimplexId, 3> triangleVec
            = {vertexIds[1], vertexIds[2], vertexIds[3]};
          SimplexId nodeId = vertexIndices_[vertexIds[1]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            nodeMaps[nodeId] = ImplicitCluster(nodeId);
            buildInternalTriangleMap(&nodeMaps[nodeId], false, true);
          }
          SimplexId localVertexId
            = vertexIds[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
              + triangleIntervals_[nodeId - 1];
          linksCount[localVertexId]++;
          // v2: (v1, v3, v4)
          if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[0] = vertexIds[0];
            localVertexId
              = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodePtr->internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[nodePtr->nid - 1];
            linksCount[localVertexId]++;
            // v3: (v1, v2, v4)
            if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
              triangleVec[1] = vertexIds[1];
              localVertexId
                = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
              vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
                = nodePtr->internalTriangleMap_.at(triangleVec)
                  + triangleIntervals_[nodePtr->nid - 1];
              linksCount[localVertexId]++;
            }
            // v4: (v1, v2, v3)
            if(vertexIds[3] <= vertexIntervals_[nodePtr->nid]) {
              triangleVec[2] = vertexIds[2];
              localVertexId
                = vertexIds[3] - vertexIntervals_[nodePtr->nid - 1] - 1;
              vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
                = nodePtr->internalTriangleMap_.at(triangleVec)
                  + triangleIntervals_[nodePtr->nid - 1];
              linksCount[localVertexId]++;
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          std::array<SimplexId, 4> vertexIds
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2),
               (SimplexId)cellArray_->getCellVertex(cid, 3)};
          // start from v2
          std::array<SimplexId, 3> triangleVec
            = {vertexIds[0], vertexIds[2], vertexIds[3]};
          SimplexId nodeId = vertexIndices_[vertexIds[0]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            nodeMaps[nodeId] = ImplicitCluster(nodeId);
            buildInternalTriangleMap(&nodeMaps[nodeId], false, true);
          }
          SimplexId localVertexId
            = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
          if(vertexIds[1] > vertexIntervals_[nodePtr->nid - 1]
             && vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[nodeId - 1];
            linksCount[localVertexId]++;
          }
          if(vertexIds[2] > vertexIntervals_[nodePtr->nid - 1]
             && vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[1] = vertexIds[1];
            localVertexId
              = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[nodeId - 1];
            linksCount[localVertexId]++;
          }
          if(vertexIds[3] > vertexIntervals_[nodePtr->nid - 1]
             && vertexIds[3] <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[1] = vertexIds[1];
            triangleVec[2] = vertexIds[2];
            localVertexId
              = vertexIds[3] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
                + triangleIntervals_[nodeId - 1];
            linksCount[localVertexId]++;
          }
        }

        // fill FlatJaggedArray struct
        nodePtr->vertexLinks_.setData(
          std::move(vertexLinkData), std::move(offsets));
      }

      return 0;
    }

    /**
     * Get the vertex neighbors for all the vertices in a given cluster.
     */
    int getClusterVertexNeighbors(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> vertexNeighborData, offsets(localVertexNum + 1, 0);

      SimplexId v1, v2;
      std::vector<boost::unordered_set<SimplexId>> vertexNeighborSet(
        localVertexNum);

      // loop through the internal cells
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          v1 = cellArray_->getCellVertex(cid, j);
          if(v1 <= vertexIntervals_[nodePtr->nid]) {
            for(SimplexId k = j + 1; k < verticesPerCell; k++) {
              v2 = cellArray_->getCellVertex(cid, k);
              vertexNeighborSet[v1 - vertexIntervals_[nodePtr->nid - 1] - 1]
                .insert(v2);
              if(v2 <= vertexIntervals_[nodePtr->nid]) {
                vertexNeighborSet[v2 - vertexIntervals_[nodePtr->nid - 1] - 1]
                  .insert(v1);
              }
            }
          }
        }
      }

      // loop through external cells
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            v1 = cellArray_->getCellVertex(cid, j);
            v2 = cellArray_->getCellVertex(cid, k);
            if(v1 > vertexIntervals_[nodePtr->nid - 1]
               && v1 <= vertexIntervals_[nodePtr->nid])
              vertexNeighborSet[v1 - vertexIntervals_[nodePtr->nid - 1] - 1]
                .insert(v2);
            if(v2 > vertexIntervals_[nodePtr->nid - 1]
               && v2 <= vertexIntervals_[nodePtr->nid])
              vertexNeighborSet[v2 - vertexIntervals_[nodePtr->nid - 1] - 1]
                .insert(v1);
          }
        }
      }

      for(SimplexId i = 1; i <= localVertexNum; i++) {
        offsets[i] = offsets[i - 1] + vertexNeighborSet[i - 1].size();
        vertexNeighborData.insert(vertexNeighborData.end(),
                                  vertexNeighborSet[i - 1].begin(),
                                  vertexNeighborSet[i - 1].end());
      }

      nodePtr->vertexNeighbors_.setData(
        std::move(vertexNeighborData), std::move(offsets));
      return 0;
    }

    /**
     * Get the vertex stars for all the vertices in a given cluster.
     * The function is similar as getVertexCells().
     */
    int getClusterVertexStars(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localVertexNum + 1, 0),
        starsCount(localVertexNum, 0);

      // set the offsets vector
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        SimplexId vertexId = cellArray_->getCellVertex(cid, 0);
        offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
        for(SimplexId j = 1; j < verticesPerCell; j++) {
          vertexId = cellArray_->getCellVertex(cid, j);
          if(vertexId > vertexIntervals_[nodePtr->nid - 1]
             && vertexId <= vertexIntervals_[nodePtr->nid]) {
            offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
          }
        }
      }
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        for(SimplexId j = 1; j < verticesPerCell; j++) {
          SimplexId vertexId = cellArray_->getCellVertex(cid, j);
          if(vertexId > vertexIntervals_[nodePtr->nid - 1]
             && vertexId <= vertexIntervals_[nodePtr->nid]) {
            offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
          }
        }
      }

      // compute partial sum of number of stars per vertex
      for(SimplexId i = 1; i <= localVertexNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for vertex star data
      std::vector<SimplexId> vertexStarData(offsets.back());

      // fill the flat vector using offsets and count vectors
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        SimplexId localVertexId = cellArray_->getCellVertex(cid, 0)
                                  - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
          = cid;
        starsCount[localVertexId]++;
        for(SimplexId j = 1; j < verticesPerCell; j++) {
          SimplexId vertexId = cellArray_->getCellVertex(cid, j);
          // see if it is in the current node
          if(vertexId > vertexIntervals_[nodePtr->nid - 1]
             && vertexId <= vertexIntervals_[nodePtr->nid]) {
            localVertexId = vertexId - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
              = cid;
            starsCount[localVertexId]++;
          }
        }
      }
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        for(SimplexId j = 0; j < verticesPerCell; j++) {
          // see if it is in the current node
          SimplexId vertexId = cellArray_->getCellVertex(cid, j);
          if(vertexId > vertexIntervals_[nodePtr->nid - 1]
             && vertexId <= vertexIntervals_[nodePtr->nid]) {
            SimplexId localVertexId
              = vertexId - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
              = cid;
            starsCount[localVertexId]++;
          }
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->vertexStars_.setData(
        std::move(vertexStarData), std::move(offsets));

      // sort the vector for cell neighbor relation
      // for(size_t i = 0; i < localVertexStars.size(); i++) {
      //   sort(localVertexStars[i].begin(), localVertexStars[i].end());
      // }

      return 0;
    }

    /**
     * Get the vertex triangles for all the vertices in a given cluster.
     */
    int getClusterVertexTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      std::vector<SimplexId> offsets(localVertexNum + 1, 0),
        trianglesCount(localVertexNum, 0);

      if(nodePtr->internalTriangleMap_.empty()) {
        buildInternalTriangleMap(nodePtr, false, true);
      }
      if(nodePtr->externalTriangleMap_.empty()) {
        buildExternalTriangleMap(nodePtr);
      }

      // set the offsets vector
      boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[j] <= vertexIntervals_[nodePtr->nid])
            offsets[iter->first[j] - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
      for(iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[j] <= vertexIntervals_[nodePtr->nid])
            offsets[iter->first[j] - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }

      // compute partial sum of number of triangles per vertex
      for(SimplexId i = 1; i <= localVertexNum; i++) {
        offsets[i] += offsets[i - 1];
      }

      // allocate the flat vector for vertex triangle data
      std::vector<SimplexId> vertexTriangleData(offsets.back());

      // fill the flat vector using offsets and count vectors
      for(iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[j] <= vertexIntervals_[nodePtr->nid]) {
            SimplexId localVertexId
              = iter->first[j] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexTriangleData[offsets[localVertexId]
                               + trianglesCount[localVertexId]]
              = iter->second + triangleIntervals_[nodePtr->nid - 1];
            trianglesCount[localVertexId]++;
          }
        }
      }
      for(iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[j] <= vertexIntervals_[nodePtr->nid]) {
            SimplexId localVertexId
              = iter->first[j] - vertexIntervals_[nodePtr->nid - 1] - 1;
            vertexTriangleData[offsets[localVertexId]
                               + trianglesCount[localVertexId]]
              = iter->second;
            trianglesCount[localVertexId]++;
          }
        }
      }

      // fill FlatJaggedArray struct
      nodePtr->vertexTriangles_.setData(
        std::move(vertexTriangleData), std::move(offsets));

      return 0;
    }

    /**
     * Get the boundary cells in a given cluster.
     */
    int getBoundaryCells(ImplicitCluster *const nodePtr,
                         const SimplexId dim = 2) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      if(getDimensionality() == 2) {
        SimplexId localEdgeNum
          = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
        if(nodePtr->boundaryEdges_.empty()) {
          nodePtr->boundaryEdges_ = std::vector<bool>(localEdgeNum, false);
          if(nodePtr->edgeStars_.empty()) {
            getClusterEdgeStars(nodePtr);
          }
          for(SimplexId i = 0; i < localEdgeNum; i++) {
            if(nodePtr->edgeStars_.size(i) == 1) {
              nodePtr->boundaryEdges_.at(i) = true;
            }
          }
        }
        // if boundary vertices are requested
        if(dim == 0 && nodePtr->boundaryVertices_.empty()) {
          nodePtr->boundaryVertices_ = std::vector<bool>(
            vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalEdgeMap_.empty()) {
            buildExternalEdgeMap(nodePtr);
          }
          // internal edges
          for(auto iter = nodePtr->internalEdgeMap_.begin();
              iter != nodePtr->internalEdgeMap_.end(); iter++) {
            if((nodePtr->boundaryEdges_)[iter->second - 1]) {
              (nodePtr
                 ->boundaryVertices_)[iter->first[0]
                                      - vertexIntervals_[nodePtr->nid - 1] - 1]
                = true;
              if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                (nodePtr->boundaryVertices_)
                  [iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
            }
          }
          // external edges
          boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
          for(auto iter = nodePtr->externalEdgeMap_.begin();
              iter != nodePtr->externalEdgeMap_.end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = ImplicitCluster(nodeId);
              getBoundaryCells(&nodeMaps[nodeId]);
            }
            if((nodeMaps[nodeId]
                  .boundaryEdges_)[iter->second - edgeIntervals_[nodeId - 1]
                                   - 1]) {
              (nodePtr
                 ->boundaryVertices_)[iter->first[1]
                                      - vertexIntervals_[nodePtr->nid - 1] - 1]
                = true;
            }
          }
        }
      } else if(getDimensionality() == 3) {
        // get the boundary triangles first
        SimplexId localTriangleNum = triangleIntervals_[nodePtr->nid]
                                     - triangleIntervals_[nodePtr->nid - 1];
        if(nodePtr->boundaryTriangles_.empty()) {
          nodePtr->boundaryTriangles_
            = std::vector<bool>(localTriangleNum, false);
          if(nodePtr->triangleStars_.empty()) {
            getClusterTriangleStars(nodePtr);
          }
          for(SimplexId i = 0; i < localTriangleNum; i++) {
            if(nodePtr->triangleStars_.size(i) == 1) {
              (nodePtr->boundaryTriangles_)[i] = true;
            }
          }
        }
        // if the boundary edges are requested
        if(dim == 1 && nodePtr->boundaryEdges_.empty()) {
          nodePtr->boundaryEdges_ = std::vector<bool>(
            edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalTriangleMap_.empty()) {
            buildExternalTriangleMap(nodePtr);
          }
          if(nodePtr->internalEdgeMap_.empty()) {
            buildInternalEdgeMap(nodePtr, false, true);
          }
          // internal triangles
          for(auto iter = nodePtr->internalTriangleMap_.begin();
              iter != nodePtr->internalTriangleMap_.end(); iter++) {
            if((nodePtr->boundaryTriangles_)[iter->second - 1]) {
              std::array<SimplexId, 2> edgePair
                = {iter->first[0], iter->first[1]};
              (nodePtr
                 ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
                = true;
              edgePair[1] = iter->first[2];
              (nodePtr
                 ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
                = true;
              if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                edgePair[0] = iter->first[1];
                (nodePtr
                   ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
                  = true;
              }
            }
          }
          // external triangles
          boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
          for(auto iter = nodePtr->externalTriangleMap_.begin();
              iter != nodePtr->externalTriangleMap_.end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = ImplicitCluster(nodeId);
              getBoundaryCells(&nodeMaps[nodeId]);
            }
            if((nodeMaps[nodeId]
                  .boundaryTriangles_)[iter->second
                                       - triangleIntervals_[nodeId - 1] - 1]) {
              if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                std::array<SimplexId, 2> edgePair
                  = {iter->first[1], iter->first[2]};
                (nodePtr
                   ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
                  = true;
              }
            }
          }
        }

        // if the boundary vertices are requested
        else if(dim == 0 && nodePtr->boundaryVertices_.empty()) {
          nodePtr->boundaryVertices_ = std::vector<bool>(
            vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalTriangleMap_.empty()) {
            buildExternalTriangleMap(nodePtr);
          }
          // internal triangles
          for(auto iter = nodePtr->internalTriangleMap_.begin();
              iter != nodePtr->internalTriangleMap_.end(); iter++) {
            if((nodePtr->boundaryTriangles_)[iter->second - 1]) {
              for(int j = 0; j < 3; j++) {
                SimplexId vid = iter->first[j];
                if(vid <= vertexIntervals_[nodePtr->nid]) {
                  (nodePtr->boundaryVertices_)
                    [vid - vertexIntervals_[nodePtr->nid - 1] - 1]
                    = true;
                }
              }
            }
          }
          // external triangles
          boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
          for(auto iter = nodePtr->externalTriangleMap_.begin();
              iter != nodePtr->externalTriangleMap_.end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = ImplicitCluster(nodeId);
              ;
              getBoundaryCells(&nodeMaps[nodeId]);
            }
            if((nodeMaps[nodeId]
                  .boundaryTriangles_)[iter->second
                                       - triangleIntervals_[nodeId - 1] - 1]) {
              if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                (nodePtr->boundaryVertices_)
                  [iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
              if(iter->first[2] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[2] <= vertexIntervals_[nodePtr->nid]) {
                (nodePtr->boundaryVertices_)
                  [iter->first[2] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
            }
          }
        }
      } else {
        return -1;
      }

      return 0;
    }

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