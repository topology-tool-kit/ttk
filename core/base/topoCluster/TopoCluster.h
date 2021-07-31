/// \ingroup base
/// \class ttk::TopoCluster
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date November 2021.
///
/// \brief TTK processing package for the TopoCluster data structure.
///
/// %TopoCluster is a localized data structure for simplicial meshes.
/// The key idea of TopoCluster is to subdivide the simplicial mesh into
/// clusters. Then, the connectivity information is computed locally for each
/// cluster and discarded when it is no longer needed. The simplicial mesh needs
/// to be subdivided before using TopoCluster, i.e., there has to be a scalar
/// field named "_index" to denote the cluster index of each vertex. Note
/// Topocluster will reindex the simplices based on the clustering input array.
///
/// \sa ttk::Triangulation
/// \sa ttk::PreTopoCluster

#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CellArray.h>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <list>

#define EDGE_ID 1
#define TRIANGLE_ID 2

using namespace std;
namespace ttk {

  class ImplicitCluster {
  private:
    /* components */
    SimplexId nid;
    vector<array<SimplexId, 2>> *internalEdgeList_;
    vector<array<SimplexId, 3>> *internalTriangleList_;
    boost::unordered_map<array<SimplexId, 2>, SimplexId> *internalEdgeMap_;
    boost::unordered_map<array<SimplexId, 2>, SimplexId> *externalEdgeMap_;
    boost::unordered_map<array<SimplexId, 3>, SimplexId> *internalTriangleMap_;
    boost::unordered_map<array<SimplexId, 3>, SimplexId> *externalTriangleMap_;
    /* boundary cells */
    vector<bool> *boundaryVertices_;
    vector<bool> *boundaryEdges_;
    vector<bool> *boundaryTriangles_;
    /* vertex relationships */
    vector<vector<SimplexId>> *vertexEdges_;
    vector<vector<SimplexId>> *vertexLinks_;
    vector<vector<SimplexId>> *vertexNeighbors_;
    vector<vector<SimplexId>> *vertexStars_;
    vector<vector<SimplexId>> *vertexTriangles_;
    /* edge relationships */
    // edgeVertex relation can be extracted from internal edge list
    vector<vector<SimplexId>> *edgeLinks_;
    vector<vector<SimplexId>> *edgeStars_;
    vector<vector<SimplexId>> *edgeTriangles_;
    /* triangle relationships */
    // triangleVertex relation can be extracted from internal triangle list
    vector<vector<SimplexId>> *triangleEdges_;
    vector<vector<SimplexId>> *triangleLinks_;
    vector<vector<SimplexId>> *triangleStars_;
    /* cell relationships */
    vector<vector<SimplexId>> *cellEdges_;
    vector<vector<SimplexId>> *cellNeighbors_;
    vector<vector<SimplexId>> *cellTriangles_;

  public:
    ImplicitCluster(SimplexId id) {
      /* components */
      nid = id;
      internalEdgeList_ = nullptr;
      internalTriangleList_ = nullptr;
      internalEdgeMap_ = nullptr;
      externalEdgeMap_ = nullptr;
      internalTriangleMap_ = nullptr;
      externalTriangleMap_ = nullptr;
      /* boundary cells */
      boundaryEdges_ = nullptr;
      boundaryVertices_ = nullptr;
      boundaryTriangles_ = nullptr;
      /* vertex relationships */
      vertexEdges_ = nullptr;
      vertexLinks_ = nullptr;
      vertexNeighbors_ = nullptr;
      vertexStars_ = nullptr;
      vertexTriangles_ = nullptr;
      /* edge relationships */
      edgeLinks_ = nullptr;
      edgeStars_ = nullptr;
      edgeTriangles_ = nullptr;
      /* triangle relationships */
      triangleLinks_ = nullptr;
      triangleEdges_ = nullptr;
      triangleStars_ = nullptr;
      /* cell relationships */
      cellEdges_ = nullptr;
      cellNeighbors_ = nullptr;
      cellTriangles_ = nullptr;
    }
    ~ImplicitCluster() {
      delete internalEdgeList_;
      delete internalTriangleList_;
      delete internalEdgeMap_;
      delete externalEdgeMap_;
      delete internalTriangleMap_;
      delete externalTriangleMap_;
      delete boundaryEdges_;
      delete boundaryVertices_;
      delete boundaryTriangles_;
      delete vertexEdges_;
      delete vertexLinks_;
      delete vertexNeighbors_;
      delete vertexStars_;
      delete vertexTriangles_;
      delete edgeLinks_;
      delete edgeStars_;
      delete edgeTriangles_;
      delete triangleEdges_;
      delete triangleLinks_;
      delete triangleStars_;
      delete cellEdges_;
      delete cellNeighbors_;
      delete cellTriangles_;
    }

    friend class TopoCluster;
  };

  class TopoCluster final : public AbstractTriangulation {

  public:
    TopoCluster();

    ~TopoCluster();

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
      vector<SimplexId> vertexMap(vertexNumber_);
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
        vector<SimplexId> vertexMap(vertexNumber_);
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
    int reorderVertices(vector<SimplexId> &vertexMap) {
      // get the number of nodes (the max value in the array)
      for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
        if(vertexIndices_[vid] > nodeNumber_) {
          nodeNumber_ = vertexIndices_[vid];
        }
      }
      nodeNumber_++; // since the index starts from 0
      vector<vector<SimplexId>> nodeVertices(nodeNumber_);
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
        double *newPointSet = new double[3 * vertexNumber_];
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
        delete[] newPointSet;
      } else {
        float *newPointSet = new float[3 * vertexNumber_];
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
        delete[] newPointSet;
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
    int reorderCells(const vector<SimplexId> &vertexMap,
                     const SimplexId &cellNumber,
                     const LongSimplexId *connectivity,
                     const LongSimplexId *offset) {
      // change the indices in cell array
      SimplexId cellCount = 0, verticesPerCell = offset[1] - offset[0];
      vector<vector<SimplexId>> nodeCells(nodeNumber_ + 1);
      LongSimplexId *cellArr = ((LongSimplexId *)connectivity);

      for(SimplexId cid = 0; cid < cellNumber; cid++) {
        SimplexId cellId = offset[cid];
        for(int j = 0; j < verticesPerCell; j++) {
          cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
        }
        sort(cellArr + cellId, cellArr + cellId + verticesPerCell);
        nodeCells[vertexIndices_[cellArr[cellId]]].push_back(cid);
      }

      // rearange the cell array
      cellIntervals_.resize(nodeNumber_ + 1);
      externalCells_.resize(nodeNumber_ + 1);
      cellIntervals_[0] = -1;
      LongSimplexId *newCellArray
        = new LongSimplexId[verticesPerCell * cellNumber];
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
      delete[] newCellArray;

      return 0;
    }

    /**
     * Reorder the input cells.
     */
    int reorderCells(const vector<SimplexId> &vertexMap,
                     const LongSimplexId *cellArray) {
      // change the indices in cell array
      SimplexId cellCount = 0, verticesPerCell = cellArray[0];
      vector<vector<SimplexId>> nodeCells(nodeNumber_ + 1);
      LongSimplexId *cellArr = ((LongSimplexId *)cellArray);

      for(SimplexId cid = 0; cid < cellNumber_; cid++) {
        SimplexId cellId = (verticesPerCell + 1) * cid;
        for(int j = 1; j <= verticesPerCell; j++) {
          cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
        }
        sort(cellArr + cellId + 1, cellArr + cellId + 1 + verticesPerCell);
        nodeCells[vertexIndices_[cellArr[cellId + 1]]].push_back(cid);
      }

      // rearange the cell array
      cellIntervals_.resize(nodeNumber_ + 1);
      externalCells_.resize(nodeNumber_ + 1);
      cellIntervals_[0] = -1;
      LongSimplexId *newCellArray
        = new LongSimplexId[(verticesPerCell + 1) * cellNumber_];
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
      delete[] newCellArray;

      return 0;
    }

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localEdgeId < 0))
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->cellEdges_ == nullptr) {
        exnode->cellEdges_ = new vector<vector<SimplexId>>();
        getCellEdges(exnode);
      }
      if(localEdgeId >= (SimplexId)(*(exnode->cellEdges_))[localCellId].size())
        return -2;

      edgeId = (*(exnode->cellEdges_))[localCellId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      return (maxCellDim_ + 1) * maxCellDim_ / 2;
    }

    inline const vector<vector<SimplexId>> *getCellEdgesInternal() override {
      cellEdgeVector_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->cellEdges_ == nullptr) {
          exnode->cellEdges_ = new vector<vector<SimplexId>>();
          getCellEdges(exnode);
        }
        cellEdgeVector_.insert(cellEdgeVector_.end(),
                               exnode->cellEdges_->begin(),
                               exnode->cellEdges_->end());
      }
      return &cellEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localNeighborId < 0))
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->cellNeighbors_ == nullptr) {
        exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
        getCellNeighbors(exnode);
      }
      if(localNeighborId
         >= (SimplexId)(*(exnode->cellNeighbors_))[localCellId].size())
        return -2;

      neighborId = (*(exnode->cellNeighbors_))[localCellId][localNeighborId];
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
      if(exnode->cellNeighbors_ == nullptr) {
        exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
        getCellNeighbors(exnode);
      }
      return (*(exnode->cellNeighbors_))[localCellId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborList_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->cellNeighbors_ == nullptr) {
          exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
          getCellNeighbors(exnode);
        }
        cellNeighborList_.insert(cellNeighborList_.end(),
                                 exnode->cellNeighbors_->begin(),
                                 exnode->cellNeighbors_->end());
      }
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= getCellTriangleNumber(cellId)))
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->cellTriangles_ == nullptr) {
        exnode->cellTriangles_ = new vector<vector<SimplexId>>();
        getCellTriangles(exnode);
      }
      triangleId = (*(exnode->cellTriangles_))[localCellId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      return (maxCellDim_ + 1) * maxCellDim_ * (maxCellDim_ - 1) / 6;
    }

    inline const vector<vector<SimplexId>> *
      getCellTrianglesInternal() override {
      cellTriangleVector_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->cellTriangles_ == nullptr) {
          exnode->cellTriangles_ = new vector<vector<SimplexId>>();
          getCellTriangles(exnode);
        }
        cellTriangleVector_.insert(cellTriangleVector_.end(),
                                   exnode->cellTriangles_->begin(),
                                   exnode->cellTriangles_->end());
      }
      return &cellTriangleVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localVertexId < 0)
         || (localVertexId >= cellArray_->getCellVertexNumber(cellId)))
        return -2;
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

    inline const vector<array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      edgeList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        vector<array<SimplexId, 2>> localInternalEdges;
        boost::unordered_map<array<SimplexId, 2>, SimplexId>
          localInternalEdgeMap;
        buildInternalEdgeMap(
          nid, &localInternalEdges, &localInternalEdgeMap, nullptr);
        edgeList_.insert(edgeList_.end(), localInternalEdges.begin(),
                         localInternalEdges.end());
      }
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeLinks_ == nullptr) {
        exnode->edgeLinks_ = new vector<vector<SimplexId>>();
        getEdgeLinks(exnode);
      }
      if(localLinkId >= (SimplexId)(*(exnode->edgeLinks_))[localEdgeId].size())
        return -2;

      linkId = (*(exnode->edgeLinks_))[localEdgeId][localLinkId];
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
      if(exnode->edgeLinks_ == nullptr) {
        exnode->edgeLinks_ = new vector<vector<SimplexId>>();
        getEdgeLinks(exnode);
      }
      return (*(exnode->edgeLinks_))[localEdgeId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeLinks_ == nullptr) {
          exnode->edgeLinks_ = new vector<vector<SimplexId>>();
          getEdgeLinks(exnode);
        }
        edgeLinkList_.insert(edgeLinkList_.end(), exnode->edgeLinks_->begin(),
                             exnode->edgeLinks_->end());
      }
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeStars_ == nullptr) {
        exnode->edgeStars_ = new vector<vector<SimplexId>>();
        buildInternalEdgeMap(nid, nullptr, nullptr, exnode->edgeStars_);
      }
      if(localStarId >= (SimplexId)(*(exnode->edgeStars_))[localEdgeId].size())
        return -2;
      starId = (*(exnode->edgeStars_))[localEdgeId][localStarId];
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
      if(exnode->edgeStars_ == nullptr) {
        exnode->edgeStars_ = new vector<vector<SimplexId>>();
        buildInternalEdgeMap(nid, nullptr, nullptr, exnode->edgeStars_);
      }
      return (*(exnode->edgeStars_))[localEdgeId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeStars_ == nullptr) {
          exnode->edgeStars_ = new vector<vector<SimplexId>>();
          buildInternalEdgeMap(nid, nullptr, nullptr, exnode->edgeStars_);
        }
        edgeStarList_.insert(edgeStarList_.end(), exnode->edgeStars_->begin(),
                             exnode->edgeStars_->end());
      }
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
      if(localTriangleId < 0)
        return -2;
#endif

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeTriangles_ == nullptr) {
        exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
        getEdgeTriangles(exnode);
      }
      if(localTriangleId
         >= (SimplexId)(*(exnode->edgeTriangles_))[localEdgeId].size())
        return -2;
      triangleId = (*(exnode->edgeTriangles_))[localEdgeId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->edgeTriangles_ == nullptr) {
        exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
        getEdgeTriangles(exnode);
      }
      return (*(exnode->edgeTriangles_))[edgeId - edgeIntervals_[nid - 1] - 1]
        .size();
    }

    inline const vector<vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->edgeTriangles_ == nullptr) {
          exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
          getEdgeTriangles(exnode);
        }
        edgeTriangleList_.insert(edgeTriangleList_.end(),
                                 exnode->edgeTriangles_->begin(),
                                 exnode->edgeTriangles_->end());
      }
      return &edgeTriangleList_;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
      if((localVertexId != 0) && (localVertexId != 1))
        return -2;
#endif

      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->internalEdgeList_ == nullptr) {
        exnode->internalEdgeList_ = new vector<array<SimplexId, 2>>();
        if(exnode->internalEdgeMap_ == nullptr) {
          exnode->internalEdgeMap_
            = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
        }
        buildInternalEdgeMap(
          nid, exnode->internalEdgeList_, exnode->internalEdgeMap_, nullptr);
      }

      if(localVertexId) {
        vertexId = (*(exnode->internalEdgeList_))[localEdgeId][1];
      } else {
        vertexId = (*(exnode->internalEdgeList_))[localEdgeId][0];
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

    inline const vector<array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override {
      // if it is a triangle mesh
      if(getDimensionality() == 2) {
        triangleList_.resize(cellNumber_, array<SimplexId, 3>());
        for(SimplexId cid = 0; cid < cellNumber_; cid++) {
          triangleList_[cid][0] = cellArray_->getCellVertex(cid, 0);
          triangleList_[cid][1] = cellArray_->getCellVertex(cid, 1);
          triangleList_[cid][2] = cellArray_->getCellVertex(cid, 2);
        }
      } else {
        triangleList_.reserve(triangleIntervals_.back() + 1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          ImplicitCluster *exnode = searchCache(nid);
          if(exnode->internalTriangleList_ == nullptr) {
            exnode->internalTriangleList_ = new vector<array<SimplexId, 3>>();
            if(exnode->internalTriangleMap_ == nullptr) {
              exnode->internalTriangleMap_
                = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            }
            buildInternalTriangleMap(nid, exnode->internalTriangleList_,
                                     exnode->internalTriangleMap_, nullptr);
          }
          triangleList_.insert(triangleList_.end(),
                               exnode->internalTriangleList_->begin(),
                               exnode->internalTriangleList_->end());
        }
      }
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if((localEdgeId < 0) || (localEdgeId > 2))
        return -2;
#endif

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleEdges_ == nullptr) {
        exnode->triangleEdges_ = new vector<vector<SimplexId>>();
        getTriangleEdges(exnode);
      }
      edgeId = (*(exnode->triangleEdges_))[localTriangleId][localEdgeId];
      return 0;
    }

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      return 3;
    }

    inline const vector<vector<SimplexId>> *
      getTriangleEdgesInternal() override {
      triangleEdgeVector_.reserve(triangleIntervals_.size() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->triangleEdges_ == nullptr) {
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(exnode);
        }
        triangleEdgeVector_.insert(triangleEdgeVector_.end(),
                                   exnode->triangleEdges_->begin(),
                                   exnode->triangleEdges_->end());
      }
      return &triangleEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleLinks_ == nullptr) {
        exnode->triangleLinks_ = new vector<vector<SimplexId>>();
        getTriangleLinks(exnode);
      }
      if(localLinkId
         >= (SimplexId)(*(exnode->triangleLinks_))[localTriangleId].size())
        return -2;

      linkId = (*(exnode->triangleLinks_))[localTriangleId][localLinkId];
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
      if(exnode->triangleLinks_ == nullptr) {
        exnode->triangleLinks_ = new vector<vector<SimplexId>>();
        getTriangleLinks(exnode);
      }
      return (*(exnode->triangleLinks_))[localTriangleId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->triangleLinks_ == nullptr) {
          exnode->triangleLinks_ = new vector<vector<SimplexId>>();
          getTriangleLinks(exnode);
        }
        triangleLinkList_.insert(triangleLinkList_.end(),
                                 exnode->triangleLinks_->begin(),
                                 exnode->triangleLinks_->end());
      }
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || !triangleIntervals_.size()
         || (triangleId > triangleIntervals_.back()))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->triangleStars_ == nullptr) {
        exnode->triangleStars_ = new vector<vector<SimplexId>>();
        buildInternalTriangleMap(nid, nullptr, nullptr, exnode->triangleStars_);
      }
      if(localStarId
         >= (SimplexId)(*(exnode->triangleStars_))[localTriangleId].size())
        return -2;

      starId = (*(exnode->triangleStars_))[localTriangleId][localStarId];
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
      if(exnode->triangleStars_ == nullptr) {
        exnode->triangleStars_ = new vector<vector<SimplexId>>();
        buildInternalTriangleMap(nid, nullptr, nullptr, exnode->triangleStars_);
      }

      return (*(exnode->triangleStars_))[localTriangleId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->triangleStars_ == nullptr) {
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          buildInternalTriangleMap(
            nid, nullptr, nullptr, exnode->triangleStars_);
        }
        triangleStarList_.insert(triangleStarList_.end(),
                                 exnode->triangleStars_->begin(),
                                 exnode->triangleStars_->end());
      }
      return &triangleStarList_;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if((localVertexId < 0) || (localVertexId > 2))
        return -2;
#endif

      SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->internalTriangleList_ == nullptr) {
        exnode->internalTriangleList_ = new vector<array<SimplexId, 3>>();
        if(exnode->internalTriangleMap_ == nullptr) {
          exnode->internalTriangleMap_
            = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        }
        buildInternalTriangleMap(nid, exnode->internalTriangleList_,
                                 exnode->internalTriangleMap_, nullptr);
      }
      vertexId
        = (*(exnode->internalTriangleList_))[localTriangleId][localVertexId];
      return 0;
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localEdgeId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexEdges_ == nullptr) {
        exnode->vertexEdges_ = new vector<vector<SimplexId>>();
        getVertexEdges(exnode);
      }
      if(localEdgeId
         >= (SimplexId)(*exnode->vertexEdges_)[localVertexId].size())
        return -2;
      edgeId = (*(exnode->vertexEdges_))[localVertexId][localEdgeId];
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
      if(exnode->vertexEdges_ == nullptr) {
        exnode->vertexEdges_ = new vector<vector<SimplexId>>();
        getVertexEdges(exnode);
      }
      return (*(exnode->vertexEdges_))[localVertexId].size();
    }

    inline const vector<vector<SimplexId>> *getVertexEdgesInternal() override {
      vertexEdgeList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexEdges_ == nullptr) {
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(exnode);
        }
        vertexEdgeList_.insert(vertexEdgeList_.end(),
                               exnode->vertexEdges_->begin(),
                               exnode->vertexEdges_->end());
      }

      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId > vertexIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexLinks_ == nullptr) {
        exnode->vertexLinks_ = new vector<vector<SimplexId>>();
        getVertexLinks(exnode);
      }
      if(localLinkId
         >= (SimplexId)(*(exnode->vertexLinks_))[localVertexId].size())
        return -2;

      linkId = (*(exnode->vertexLinks_))[localVertexId][localLinkId];
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
      if(exnode->vertexLinks_ == nullptr) {
        exnode->vertexLinks_ = new vector<vector<SimplexId>>();
        getVertexLinks(exnode);
      }
      return (*(exnode->vertexLinks_))[localVertexId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkList_.reserve(vertexIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexLinks_ == nullptr) {
          exnode->vertexLinks_ = new vector<vector<SimplexId>>();
          getVertexLinks(exnode);
        }
        vertexLinkList_.insert(vertexLinkList_.end(),
                               exnode->vertexLinks_->begin(),
                               exnode->vertexLinks_->end());
      }
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localNeighborId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexNeighbors_ == nullptr) {
        exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
        getVertexNeighbors(exnode);
      }
      if(localNeighborId
         >= (SimplexId)(*(exnode->vertexNeighbors_))[localVertexId].size())
        return -2;
      neighborId
        = (*(exnode->vertexNeighbors_))[localVertexId][localNeighborId];
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
      if(exnode->vertexNeighbors_ == nullptr) {
        exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
        getVertexNeighbors(exnode);
      }
      return (*(exnode->vertexNeighbors_))[localVertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexNeighbors_ == nullptr) {
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(exnode);
        }
        vertexNeighborList_.insert(vertexNeighborList_.end(),
                                   exnode->vertexNeighbors_->begin(),
                                   exnode->vertexNeighbors_->end());
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
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexStars_ == nullptr) {
        exnode->vertexStars_ = new vector<vector<SimplexId>>();
        getVertexStars(exnode);
      }
      if(localStarId
         >= (SimplexId)(*exnode->vertexStars_)[localVertexId].size())
        return -2;
      starId = (*exnode->vertexStars_)[localVertexId][localStarId];
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
      if(exnode->vertexStars_ == nullptr) {
        exnode->vertexStars_ = new vector<vector<SimplexId>>();
        getVertexStars(exnode);
      }
      return (*(exnode->vertexStars_))[localVertexId].size();
    }

    inline const vector<vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexStars_ == nullptr) {
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(exnode);
        }
        vertexStarList_.insert(vertexStarList_.end(),
                               exnode->vertexStars_->begin(),
                               exnode->vertexStars_->end());
      }
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localTriangleId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexTriangles_ == nullptr) {
        exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
        getVertexTriangles(exnode);
      }
      if(localTriangleId
         >= (SimplexId)(*(exnode->vertexTriangles_))[localVertexId].size())
        return -2;
      triangleId
        = (*(exnode->vertexTriangles_))[localVertexId][localTriangleId];

      return 0;
    }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      ImplicitCluster *exnode = searchCache(nid);
      if(exnode->vertexTriangles_ == nullptr) {
        exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
        getVertexTriangles(exnode);
      }
      return (*(
        exnode->vertexTriangles_))[vertexId - vertexIntervals_[nid - 1] - 1]
        .size();
    }

    inline const vector<vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        ImplicitCluster *exnode = searchCache(nid);
        if(exnode->vertexTriangles_ == nullptr) {
          exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
          getVertexTriangles(exnode);
        }
        vertexTriangleList_.insert(vertexTriangleList_.end(),
                                   exnode->vertexTriangles_->begin(),
                                   exnode->vertexTriangles_->end());
      }
      return &vertexTriangleList_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return false;
#endif
      SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
      SimplexId localedgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ImplicitCluster *exnode = searchCache(nid);
      getBoundaryCells(exnode, 1);
      return (*(exnode->boundaryEdges_))[localedgeId];
    }

    bool isEmpty() const {
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
      return (*(exnode->boundaryTriangles_))[localtriangleId];
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
      return (*(exnode->boundaryVertices_))[localVertexId];
    }

    inline int preconditionBoundaryEdgesInternal() override {
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionEdgesInternal();
        hasPreconditionedBoundaryEdges_ = true;
      } else {
        // // unsupported dimension
        // std::stringstream msg;
        // msg << "[TopoCluster] Unsupported dimension for boundary "
        //   << "precondition." << std::endl;
        // dMsg(std::cerr, msg.str(), infoMsg);
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionTrianglesInternal();
        hasPreconditionedBoundaryTriangles_ = true;
      } else {
        // unsupported dimension
        // std::stringstream msg;
        // msg << "[TopoCluster] Unsupported dimension for boundary "
        //   << "precondition." << std::endl;
        // dMsg(std::cerr, msg.str(), infoMsg);
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
        vector<SimplexId> edgeCount(nodeNumber_ + 1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeCount[nid] = countInternalEdges(nid);
        }

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeIntervals_[nid] = edgeIntervals_[nid - 1] + edgeCount[nid];
        }

        cout << "[TopoCluster] Edges preconditioned in " << t.getElapsedTime()
             << " s.\n";
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() override {
      if(getDimensionality() == 2 || getDimensionality() == 3) {
        preconditionEdges();
        hasPreconditionedEdgeLinks_ = true;
      } else {
        // unsupported dimension
        // std::stringstream msg;
        // msg
        //   << "[TopoCluster] Unsupported dimension for edge link "
        //   << "precondition." << std::endl;
        // dMsg(std::cerr, msg.str(), infoMsg);
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
        std::cout << "[TopoCluster] Triangles preconditioned in "
                  << t.getElapsedTime() << " s.\n";
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
        // unsupported dimension
        // std::stringstream msg;
        // msg << "[TopoCluster] Unsupported dimension for vertex"
        //   << " link precondition." << std::endl;
        // dMsg(std::cerr, msg.str(), infoMsg);
        // return -1;
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
     * Initialize the cache with the size.
     */
    void initCache(const size_t size = 100) {
      cacheSize_ = size;
      for(int i = 0; i < threadNumber_; i++) {
        caches_[i].clear();
        cacheMaps_[i].clear();
      }
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
      const vector<SimplexId> *intervals = nullptr;
      // determine which vector to search
      if(idType == EDGE_ID) {
        intervals = &edgeIntervals_;
      } else if(idType == TRIANGLE_ID) {
        intervals = &triangleIntervals_;
      } else {
        return -1;
      }

      vector<SimplexId>::const_iterator low
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
        // missCount_++;
        if(caches_[threadId].size() >= cacheSize_) {
          if(caches_[threadId].back()->nid == reservedId) {
            return nullptr;
          }
          cacheMaps_[threadId].erase(caches_[threadId].back()->nid);
          delete caches_[threadId].back();
          caches_[threadId].pop_back();
        }
        caches_[threadId].push_front(new ImplicitCluster(nodeId));
        cacheMaps_[threadId][nodeId] = caches_[threadId].begin();
      }
      return (*cacheMaps_[threadId][nodeId]);
    }

    /**
     * Build the internal edge list in the node.
     */
    int
      buildInternalEdgeMap(SimplexId nodeId,
                           vector<array<SimplexId, 2>> *const internalEdgeList,
                           boost::unordered_map<array<SimplexId, 2>, SimplexId>
                             *const internalEdgeMap,
                           vector<vector<SimplexId>> *const edgeStars) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodeId <= 0 || nodeId > nodeNumber_)
        return -1;
#endif

      SimplexId edgeCount = 0,
                verticesPerCell = cellArray_->getCellVertexNumber(0);
      boost::unordered_map<array<SimplexId, 2>, SimplexId> *edgeMap;
      // loop through all the internal cells first
      if(edgeStars) {
        edgeStars->clear();
        edgeStars->reserve(
          6 * (vertexIntervals_[nodeId] - vertexIntervals_[nodeId - 1]));
      }

      if(internalEdgeMap) {
        // if the edge map has been computed before and only request the edge
        // list
        if(!internalEdgeMap->empty() && internalEdgeList) {
          internalEdgeList->resize(internalEdgeMap->size());
          for(auto iter = internalEdgeMap->begin();
              iter != internalEdgeMap->end(); iter++) {
            (*internalEdgeList)[iter->second - 1] = iter->first;
          }
          return 0;
        }
        internalEdgeMap->clear();
        edgeMap = internalEdgeMap;
      } else {
        edgeMap = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
      }

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        array<SimplexId, 2> edgeIds;

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
            if(edgeMap->find(edgeIds) == edgeMap->end()) {
              edgeCount++;
              (*edgeMap)[edgeIds] = edgeCount;
              if(edgeStars)
                edgeStars->push_back(vector<SimplexId>{cid});
            } else if(edgeStars) {
              edgeStars->at(edgeMap->at(edgeIds) - 1).push_back(cid);
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodeId]) {
        array<SimplexId, 2> edgeIds;

        // loop through each edge of the cell
        for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            edgeIds[0] = cellArray_->getCellVertex(cid, j);
            edgeIds[1] = cellArray_->getCellVertex(cid, k);

            // the edge is in the current node
            if(edgeIds[0] > vertexIntervals_[nodeId - 1]
               && edgeIds[0] <= vertexIntervals_[nodeId]) {
              if(edgeMap->find(edgeIds) == edgeMap->end()) {
                edgeCount++;
                (*edgeMap)[edgeIds] = edgeCount;
                if(edgeStars)
                  edgeStars->push_back(vector<SimplexId>{cid});
              } else if(edgeStars) {
                edgeStars->at(edgeMap->at(edgeIds) - 1).push_back(cid);
              }
            }
          }
        }
      }

      if(internalEdgeList) {
        internalEdgeList->resize(edgeMap->size());
        for(auto iter = edgeMap->begin(); iter != edgeMap->end(); iter++) {
          (*internalEdgeList)[iter->second - 1] = iter->first;
        }
      }

      if(!internalEdgeMap) {
        delete edgeMap;
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
      boost::unordered_map<SimplexId, vector<array<SimplexId, 2>>> edgeNodes;
      nodePtr->externalEdgeMap_
        = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();

      // loop through the external cell list
      for(size_t i = 0; i < externalCells_[nodePtr->nid].size(); i++) {
        array<SimplexId, 2> edgeIds;
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

      boost::unordered_map<SimplexId, vector<array<SimplexId, 2>>>::iterator
        iter;
      for(iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
        if(!exnode) {
          boost::unordered_map<array<SimplexId, 2>, SimplexId>
            localInternalEdgeMap;
          buildInternalEdgeMap(
            iter->first, nullptr, &localInternalEdgeMap, nullptr);
          for(array<SimplexId, 2> edgePair : iter->second) {
            (*(nodePtr->externalEdgeMap_))[edgePair]
              = localInternalEdgeMap.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          if(exnode->internalEdgeMap_ == nullptr) {
            exnode->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
          }
          for(array<SimplexId, 2> edgePair : iter->second) {
            (*(nodePtr->externalEdgeMap_))[edgePair]
              = exnode->internalEdgeMap_->at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Build the internal triangle list in the node.
     */
    int buildInternalTriangleMap(
      SimplexId nodeId,
      vector<array<SimplexId, 3>> *const internalTriangleList,
      boost::unordered_map<array<SimplexId, 3>, SimplexId>
        *const internalTriangleMap,
      vector<vector<SimplexId>> *const triangleStars) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodeId <= 0 || nodeId > nodeNumber_)
        return -1;
#endif

      SimplexId triangleCount = 0,
                verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localVertexNum
        = vertexIntervals_[nodeId] - vertexIntervals_[nodeId - 1];
      boost::unordered_map<array<SimplexId, 3>, SimplexId> *triangleMap;
      set<vector<SimplexId>> externalTriangleSet;

      if(triangleStars) {
        triangleStars->clear();
        triangleStars->reserve(9 * localVertexNum);
      }

      if(internalTriangleMap) {
        // if the triangle map has been computed before and only request the
        // triangle list
        if(!internalTriangleMap->empty() && internalTriangleList) {
          internalTriangleList->resize(internalTriangleMap->size());
          for(auto iter = internalTriangleMap->begin();
              iter != internalTriangleMap->end(); iter++) {
            (*internalTriangleList)[iter->second - 1] = iter->first;
          }
          return 0;
        }
        internalTriangleMap->clear();
        triangleMap = internalTriangleMap;
      } else {
        triangleMap
          = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
      }

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        array<SimplexId, 3> triangleIds;

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

              if(triangleMap->find(triangleIds) == triangleMap->end()) {
                triangleCount++;
                (*triangleMap)[triangleIds] = triangleCount;
                if(triangleStars) {
                  triangleStars->push_back(vector<SimplexId>{cid});
                }
              } else if(triangleStars) {
                triangleStars->at(triangleMap->at(triangleIds) - 1)
                  .push_back(cid);
              }
            }
          }
        }
      }

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodeId]) {
        array<SimplexId, 3> triangleIds;

        // loop through each triangle of the cell
        for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] > vertexIntervals_[nodeId - 1]
             && triangleIds[0] <= vertexIntervals_[nodeId]) {
            for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
              for(SimplexId l = k + 1; l < verticesPerCell; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);

                if(triangleMap->find(triangleIds) == triangleMap->end()) {
                  triangleCount++;
                  (*triangleMap)[triangleIds] = triangleCount;
                  if(triangleStars) {
                    triangleStars->push_back(vector<SimplexId>{cid});
                  }
                } else if(triangleStars) {
                  triangleStars->at(triangleMap->at(triangleIds) - 1)
                    .push_back(cid);
                }
              }
            }
          }
        }
      }

      if(internalTriangleList) {
        internalTriangleList->resize(triangleMap->size());
        for(auto iter = triangleMap->begin(); iter != triangleMap->end();
            iter++) {
          (*internalTriangleList)[iter->second - 1] = iter->first;
        }
      }

      if(!internalTriangleMap) {
        delete triangleMap;
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
      boost::unordered_map<SimplexId, vector<array<SimplexId, 3>>>
        nodeTriangles;
      nodePtr->externalTriangleMap_
        = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();

      // loop through the external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        array<SimplexId, 3> triangleIds;

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

      ImplicitCluster *exnode;
      boost::unordered_map<SimplexId, vector<array<SimplexId, 3>>>::iterator
        iter;
      for(iter = nodeTriangles.begin(); iter != nodeTriangles.end(); iter++) {
        exnode = searchCache(iter->first, nodePtr->nid);

        if(!exnode) {
          boost::unordered_map<array<SimplexId, 3>, SimplexId>
            localInternalTriangleMap;
          buildInternalTriangleMap(
            iter->first, nullptr, &localInternalTriangleMap, nullptr);
          for(array<SimplexId, 3> triangleVec : iter->second) {
            (*(nodePtr->externalTriangleMap_))[triangleVec]
              = localInternalTriangleMap.at(triangleVec)
                + triangleIntervals_[iter->first - 1];
          }
        } else {
          if(exnode->internalTriangleMap_ == nullptr) {
            exnode->internalTriangleMap_
              = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            buildInternalTriangleMap(
              iter->first, nullptr, exnode->internalTriangleMap_, nullptr);
          }
          for(array<SimplexId, 3> triangleVec : iter->second) {
            (*(nodePtr->externalTriangleMap_))[triangleVec]
              = exnode->internalTriangleMap_->at(triangleVec)
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
      boost::unordered_set<array<SimplexId, 2>,
                           boost::hash<array<SimplexId, 2>>>
        edgeSet;

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        array<SimplexId, 2> edgeIds;

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
        array<SimplexId, 2> edgeIds;

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
      boost::unordered_set<vector<SimplexId>, boost::hash<vector<SimplexId>>>
        triangleSet;

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
          cid <= cellIntervals_[nodeId]; cid++) {
        vector<SimplexId> triangleIds(3);

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
        vector<SimplexId> triangleIds(3);

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
     * Get the cell edges for all cells in a given node.
     */
    int getCellEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId edgesPerCell = verticesPerCell * (verticesPerCell - 1) / 2;
      nodePtr->cellEdges_->clear();
      nodePtr->cellEdges_->resize(
        cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1],
        vector<SimplexId>(edgesPerCell));
      boost::unordered_map<SimplexId, vector<vector<SimplexId>>> edgeNodes;

      if(nodePtr->internalEdgeMap_ == nullptr) {
        nodePtr->internalEdgeMap_
          = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
        buildInternalEdgeMap(
          nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
      }

      for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
          i <= cellIntervals_[nodePtr->nid]; i++) {
        int cnt = 0;
        // get the internal edge id from the map
        for(SimplexId k = 1; k < verticesPerCell; k++) {
          array<SimplexId, 2> edgePair
            = {(SimplexId)cellArray_->getCellVertex(i, 0),
               (SimplexId)cellArray_->getCellVertex(i, k)};
          (*(nodePtr->cellEdges_))[i - cellIntervals_[nodePtr->nid - 1] - 1]
                                  [cnt++]
            = nodePtr->internalEdgeMap_->at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
        }
        for(SimplexId j = 1; j < verticesPerCell - 1; j++) {
          for(SimplexId k = j + 1; k < verticesPerCell; k++) {
            array<SimplexId, 2> edgePair
              = {(SimplexId)cellArray_->getCellVertex(i, j),
                 (SimplexId)cellArray_->getCellVertex(i, k)};
            if(edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
              (*(nodePtr->cellEdges_))[i - cellIntervals_[nodePtr->nid - 1] - 1]
                                      [cnt++]
                = nodePtr->internalEdgeMap_->at(edgePair)
                  + edgeIntervals_[nodePtr->nid - 1];
            }
            // group the external edges by node id
            else {
              vector<SimplexId> edgeTuple{i, cnt++, edgePair[0], edgePair[1]};
              SimplexId nodeNum = vertexIndices_[edgePair[0]];
              edgeNodes[nodeNum].push_back(edgeTuple);
            }
          }
        }
      }

      ImplicitCluster *exnode;
      for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        exnode = searchCache(iter->first, nodePtr->nid);
        if(!exnode) {
          boost::unordered_map<array<SimplexId, 2>, SimplexId>
            localInternalEdgeMap;
          buildInternalEdgeMap(
            iter->first, nullptr, &localInternalEdgeMap, nullptr);
          for(vector<SimplexId> edgeTuple : iter->second) {
            array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
            (*nodePtr
                ->cellEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                              - 1][edgeTuple[1]]
              = localInternalEdgeMap.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          if(exnode->internalEdgeMap_ == nullptr) {
            exnode->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
          }
          for(vector<SimplexId> edgeTuple : iter->second) {
            array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
            (*nodePtr
                ->cellEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                              - 1][edgeTuple[1]]
              = exnode->internalEdgeMap_->at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the cell neighbors for all cells in a given node.
     */
    int getCellNeighbors(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->cellNeighbors_->clear();
      nodePtr->cellNeighbors_->resize(cellIntervals_[nodePtr->nid]
                                      - cellIntervals_[nodePtr->nid - 1]);
      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);

      if(nodePtr->vertexStars_ == nullptr) {
        nodePtr->vertexStars_ = new vector<vector<SimplexId>>();
        getVertexStars(nodePtr);
        for(int i = 0; i < (int)nodePtr->vertexStars_->size(); i++) {
          sort((*(nodePtr->vertexStars_))[i].begin(),
               (*(nodePtr->vertexStars_))[i].end());
        }
      }

      boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        for(SimplexId j = 1; j < verticesPerCell; j++) {
          if(cellArray_->getCellVertex(cid, j)
             > vertexIntervals_[nodePtr->nid]) {
            SimplexId nodeId
              = vertexIndices_[cellArray_->getCellVertex(cid, j)];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              ImplicitCluster *newNode = new ImplicitCluster(nodeId);
              newNode->vertexStars_ = new vector<vector<SimplexId>>();
              getVertexStars(newNode);
              for(int i = 0; i < (int)newNode->vertexStars_->size(); i++) {
                sort((*(newNode->vertexStars_))[i].begin(),
                     (*(newNode->vertexStars_))[i].end());
              }
              nodeMaps[nodeId] = newNode;
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

            vector<SimplexId> star0, star1;
            if(v0 <= vertexIntervals_[nodePtr->nid]) {
              star0 = (*(
                nodePtr
                  ->vertexStars_))[v0 - vertexIntervals_[nodePtr->nid - 1] - 1];
            } else {
              SimplexId nid = vertexIndices_[v0];
              star0 = (*(nodeMaps[nid]
                           ->vertexStars_))[v0 - vertexIntervals_[nid - 1] - 1];
            }
            if(v1 <= vertexIntervals_[nodePtr->nid]) {
              star1 = (*(
                nodePtr
                  ->vertexStars_))[v1 - vertexIntervals_[nodePtr->nid - 1] - 1];
            } else {
              SimplexId nid = vertexIndices_[v1];
              star1 = (*(nodeMaps[nid]
                           ->vertexStars_))[v1 - vertexIntervals_[nid - 1] - 1];
            }

            // perform an intersection of the 2 sorted star lists
            SimplexId pos0 = 0, pos1 = 0;
            SimplexId intersection = -1;

            while((pos0 < (SimplexId)star0.size())
                  && (pos1 < (SimplexId)star1.size())) {

              SimplexId biggest = star0[pos0];
              if(star1[pos1] > biggest) {
                biggest = star1[pos1];
              }

              for(SimplexId l = pos0; l < (SimplexId)star0.size(); l++) {
                if(star0[l] < biggest) {
                  pos0++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos1; l < (SimplexId)star1.size(); l++) {
                if(star1[l] < biggest) {
                  pos1++;
                } else {
                  break;
                }
              }

              if(pos0 >= (SimplexId)star0.size()
                 || pos1 >= (SimplexId)star1.size())
                break;

              if(star0[pos0] == star1[pos1]) {
                if(star0[pos0] != cid) {
                  intersection = star0[pos0];
                  break;
                }
                pos0++;
                pos1++;
              }
            }

            if(intersection != -1) {
              (*(
                nodePtr
                  ->cellNeighbors_))[cid - cellIntervals_[nodePtr->nid - 1] - 1]
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

            vector<SimplexId> star0, star1, star2;
            if(v0 <= vertexIntervals_[nodePtr->nid]) {
              star0 = (*(
                nodePtr
                  ->vertexStars_))[v0 - vertexIntervals_[nodePtr->nid - 1] - 1];
            } else {
              SimplexId nid = vertexIndices_[v0];
              star0 = (*(nodeMaps[nid]
                           ->vertexStars_))[v0 - vertexIntervals_[nid - 1] - 1];
            }
            if(v1 <= vertexIntervals_[nodePtr->nid]) {
              star1 = (*(
                nodePtr
                  ->vertexStars_))[v1 - vertexIntervals_[nodePtr->nid - 1] - 1];
            } else {
              SimplexId nid = vertexIndices_[v1];
              star1 = (*(nodeMaps[nid]
                           ->vertexStars_))[v1 - vertexIntervals_[nid - 1] - 1];
            }
            if(v2 <= vertexIntervals_[nodePtr->nid]) {
              star2 = (*(
                nodePtr
                  ->vertexStars_))[v2 - vertexIntervals_[nodePtr->nid - 1] - 1];
            } else {
              SimplexId nid = vertexIndices_[v2];
              star2 = (*(nodeMaps[nid]
                           ->vertexStars_))[v2 - vertexIntervals_[nid - 1] - 1];
            }

            // perform an intersection of the 3 (sorted) star lists
            SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
            SimplexId intersection = -1;

            while((pos0 < (SimplexId)star0.size())
                  && (pos1 < (SimplexId)star1.size())
                  && (pos2 < (SimplexId)star2.size())) {

              SimplexId biggest = star0[pos0];
              if(star1[pos1] > biggest) {
                biggest = star1[pos1];
              }
              if(star2[pos2] > biggest) {
                biggest = star2[pos2];
              }

              for(SimplexId l = pos0; l < (SimplexId)star0.size(); l++) {
                if(star0[l] < biggest) {
                  pos0++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos1; l < (SimplexId)star1.size(); l++) {
                if(star1[l] < biggest) {
                  pos1++;
                } else {
                  break;
                }
              }
              for(SimplexId l = pos2; l < (SimplexId)star2.size(); l++) {
                if(star2[l] < biggest) {
                  pos2++;
                } else {
                  break;
                }
              }

              if(pos0 >= (SimplexId)star0.size()
                 || pos1 >= (SimplexId)star1.size()
                 || pos2 >= (SimplexId)star2.size())
                break;

              if((star0[pos0] == star1[pos1]) && (star0[pos0] == star2[pos2])) {
                if(star0[pos0] != cid) {
                  intersection = star0[pos0];
                  break;
                }
                pos0++;
                pos1++;
                pos2++;
              }
            }

            if(intersection != -1) {
              (*(
                nodePtr
                  ->cellNeighbors_))[cid - cellIntervals_[nodePtr->nid - 1] - 1]
                .push_back(intersection);
            }
          }
        }
      }

      // release the memory
      for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
        delete iter->second;
      }

      return 0;
    }

    /**
     * Get the cell triangles for all cells in a given node.
     */
    int getCellTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId trianglesPerCell = getCellTriangleNumberInternal(0);
      boost::unordered_map<SimplexId, vector<vector<SimplexId>>> nodeTriangles;

      nodePtr->cellTriangles_->clear();
      nodePtr->cellTriangles_->resize(
        cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1],
        vector<SimplexId>(trianglesPerCell));

      if(nodePtr->internalTriangleMap_ == nullptr) {
        nodePtr->internalTriangleMap_
          = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        buildInternalTriangleMap(
          nodePtr->nid, nullptr, nodePtr->internalTriangleMap_, nullptr);
      }

      for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
          i <= cellIntervals_[nodePtr->nid]; i++) {
        array<SimplexId, 3> triangleVec;
        // get the internal triangle from the map
        triangleVec[0] = cellArray_->getCellVertex(i, 0);
        for(SimplexId k = 1; k < verticesPerCell - 1; k++) {
          triangleVec[1] = cellArray_->getCellVertex(i, k);
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleVec[2] = cellArray_->getCellVertex(i, l);
            (*(nodePtr->cellTriangles_))[i - cellIntervals_[nodePtr->nid - 1]
                                         - 1][k + l - 3]
              = nodePtr->internalTriangleMap_->at(triangleVec)
                + triangleIntervals_[nodePtr->nid - 1];
          }
        }
        // group the external triangles by node id
        triangleVec[0] = cellArray_->getCellVertex(i, 1);
        triangleVec[1] = cellArray_->getCellVertex(i, 2);
        triangleVec[2] = cellArray_->getCellVertex(i, 3);
        if(triangleVec[0] <= vertexIntervals_[nodePtr->nid]) {
          (*(nodePtr->cellTriangles_))[i - cellIntervals_[nodePtr->nid - 1] - 1]
            .back()
            = nodePtr->internalTriangleMap_->at(triangleVec)
              + triangleIntervals_[nodePtr->nid - 1];
        } else {
          vector<SimplexId> triangleTuple
            = {i, triangleVec[0], triangleVec[1], triangleVec[2]};
          SimplexId nodeNum = vertexIndices_[triangleVec[0]];
          nodeTriangles[nodeNum].push_back(triangleTuple);
        }
      }

      ImplicitCluster *exnode;
      for(auto iter = nodeTriangles.begin(); iter != nodeTriangles.end();
          iter++) {
        exnode = searchCache(iter->first, nodePtr->nid);
        if(!exnode) {
          boost::unordered_map<array<SimplexId, 3>, SimplexId>
            localInternalTriangleMap;
          buildInternalTriangleMap(
            iter->first, nullptr, &localInternalTriangleMap, nullptr);
          for(vector<SimplexId> triangleVec : iter->second) {
            array<SimplexId, 3> triangle
              = {triangleVec[1], triangleVec[2], triangleVec[3]};
            (*(nodePtr->cellTriangles_))[triangleVec[0]
                                         - cellIntervals_[nodePtr->nid - 1] - 1]
              .back()
              = localInternalTriangleMap.at(triangle)
                + triangleIntervals_[iter->first - 1];
          }
        } else {
          if(exnode->internalTriangleMap_ == nullptr) {
            exnode->internalTriangleMap_
              = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            buildInternalTriangleMap(
              iter->first, nullptr, exnode->internalTriangleMap_, nullptr);
          }
          for(vector<SimplexId> triangleVec : iter->second) {
            array<SimplexId, 3> triangle
              = {triangleVec[1], triangleVec[2], triangleVec[3]};
            (*(nodePtr->cellTriangles_))[triangleVec[0]
                                         - cellIntervals_[nodePtr->nid - 1] - 1]
              .back()
              = exnode->internalTriangleMap_->at(triangle)
                + triangleIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the edge links for all the edges in a given node.
     */
    int getEdgeLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->edgeLinks_->clear();
      nodePtr->edgeLinks_->resize(edgeIntervals_[nodePtr->nid]
                                  - edgeIntervals_[nodePtr->nid - 1]);
      if(getDimensionality() == 2) {
        if(nodePtr->edgeStars_ == nullptr) {
          if(nodePtr->internalEdgeMap_ == nullptr) {
            nodePtr->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
          }
          nodePtr->edgeStars_ = new vector<vector<SimplexId>>();
          buildInternalEdgeMap(nodePtr->nid, nullptr, nodePtr->internalEdgeMap_,
                               nodePtr->edgeStars_);
        }
        boost::unordered_map<array<SimplexId, 2>, SimplexId>::const_iterator
          iter;
        for(iter = nodePtr->internalEdgeMap_->begin();
            iter != nodePtr->internalEdgeMap_->end(); iter++) {
          for(SimplexId j = 0;
              j < (SimplexId)(*(nodePtr->edgeStars_))[iter->second - 1].size();
              j++) {
            SimplexId vertexId = -1;
            for(int k = 0; k < 3; k++) {
              if((cellArray_->getCellVertex(
                    (*(nodePtr->edgeStars_))[iter->second - 1][j], k)
                  != iter->first[0])
                 && (cellArray_->getCellVertex(
                       (*(nodePtr->edgeStars_))[iter->second - 1][j], k)
                     != iter->first[1])) {
                vertexId = cellArray_->getCellVertex(
                  (*(nodePtr->edgeStars_))[iter->second - 1][j], k);
                break;
              }
            }
            if(vertexId != -1) {
              (*(nodePtr->edgeLinks_))[iter->second - 1].push_back(vertexId);
            }
          }
        }
      } else if(getDimensionality() == 3) {
        if(nodePtr->cellEdges_ == nullptr) {
          nodePtr->cellEdges_ = new vector<vector<SimplexId>>();
          getCellEdges(nodePtr);
        }
        if(nodePtr->internalEdgeMap_ == nullptr) {
          nodePtr->internalEdgeMap_
            = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
          buildInternalEdgeMap(
            nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
        }

        for(SimplexId cid = 0; cid < cellIntervals_[nodePtr->nid]
                                       - cellIntervals_[nodePtr->nid - 1];
            cid++) {
          SimplexId cellId = cid + cellIntervals_[nodePtr->nid - 1] + 1;
          array<SimplexId, 2> edgePair;
          edgePair[0] = cellArray_->getCellVertex(cellId, 0);
          for(SimplexId j = 1; j < 4; j++) {
            edgePair[1] = cellArray_->getCellVertex(cellId, j);
            (*(
              nodePtr->edgeLinks_))[nodePtr->internalEdgeMap_->at(edgePair) - 1]
              .push_back((*(nodePtr->cellEdges_))[cid][6 - j]);
          }
          if(cellArray_->getCellVertex(cellId, 1)
             <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = cellArray_->getCellVertex(cellId, 1);
            for(int j = 2; j < 4; j++) {
              edgePair[1] = cellArray_->getCellVertex(cellId, j);
              (*(nodePtr
                   ->edgeLinks_))[nodePtr->internalEdgeMap_->at(edgePair) - 1]
                .push_back((*(nodePtr->cellEdges_))[cid][4 - j]);
            }
            if(cellArray_->getCellVertex(cellId, 2)
               <= vertexIntervals_[nodePtr->nid]) {
              edgePair = {(SimplexId)cellArray_->getCellVertex(cellId, 2),
                          (SimplexId)cellArray_->getCellVertex(cellId, 3)};
              (*(nodePtr
                   ->edgeLinks_))[nodePtr->internalEdgeMap_->at(edgePair) - 1]
                .push_back((*(nodePtr->cellEdges_))[cid][0]);
            }
          }
        }

        boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;
        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          array<SimplexId, 2> edgeIds;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < 3; j++) {
            for(SimplexId k = j + 1; k < 4; k++) {
              edgeIds[0] = cellArray_->getCellVertex(cid, j);
              edgeIds[1] = cellArray_->getCellVertex(cid, k);

              // the edge is in the current node
              if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
                 && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
                array<SimplexId, 2> otherEdge = {-1, -1};
                for(int i = 0; i < 4; i++) {
                  if(cellArray_->getCellVertex(cid, i) != edgeIds[0]
                     && cellArray_->getCellVertex(cid, i) != edgeIds[1]) {
                    if(otherEdge[0] == -1) {
                      otherEdge[0] = cellArray_->getCellVertex(cid, i);
                    } else if(otherEdge[1] == -1) {
                      otherEdge[1] = cellArray_->getCellVertex(cid, i);
                    } else {
                      cerr << "[TopoCluster] More than two other vertices are "
                              "found in the edge!\n";
                    }
                  }
                }
                SimplexId nodeId = vertexIndices_[otherEdge[0]];
                if(nodeMaps.find(nodeId) == nodeMaps.end()) {
                  ImplicitCluster *newNode = new ImplicitCluster(nodeId);
                  newNode->internalEdgeMap_
                    = new boost::unordered_map<array<SimplexId, 2>,
                                               SimplexId>();
                  buildInternalEdgeMap(
                    nodeId, nullptr, newNode->internalEdgeMap_, nullptr);
                  nodeMaps[nodeId] = newNode;
                }
                (*(nodePtr
                     ->edgeLinks_))[nodePtr->internalEdgeMap_->at(edgeIds) - 1]
                  .push_back(nodeMaps[nodeId]->internalEdgeMap_->at(otherEdge));
              }
            }
          }
        }

        // release the memory
        for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
          delete iter->second;
        }
      }

      return 0;
    }

    /**
     * Get the edge triangles for all the edges in a given node.
     */
    int getEdgeTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->edgeTriangles_->clear();
      nodePtr->edgeTriangles_->resize(edgeIntervals_[nodePtr->nid]
                                      - edgeIntervals_[nodePtr->nid - 1]);

      if(nodePtr->internalEdgeMap_ == nullptr) {
        nodePtr->internalEdgeMap_
          = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
        buildInternalEdgeMap(
          nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
      }
      if(nodePtr->internalTriangleMap_ == nullptr) {
        nodePtr->internalTriangleMap_
          = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        buildInternalTriangleMap(
          nodePtr->nid, nullptr, nodePtr->internalTriangleMap_, nullptr);
      }
      if(nodePtr->externalTriangleMap_ == nullptr) {
        buildExternalTriangleMap(nodePtr);
      }

      // for internal triangles
      boost::unordered_map<array<SimplexId, 3>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalTriangleMap_->begin();
          iter != nodePtr->internalTriangleMap_->end(); iter++) {
        array<SimplexId, 2> edge1 = {iter->first[0], iter->first[1]};
        array<SimplexId, 2> edge2 = {iter->first[0], iter->first[2]};

        (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge1) - 1]
          .push_back(iter->second + triangleIntervals_[nodePtr->nid - 1]);
        (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge2) - 1]
          .push_back(iter->second + triangleIntervals_[nodePtr->nid - 1]);

        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          edge1 = {iter->first[1], iter->first[2]};
          (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge1) - 1]
            .push_back(iter->second + triangleIntervals_[nodePtr->nid - 1]);
        }
      }

      // for external triangles
      // loop through each edge of the cell
      for(iter = nodePtr->externalTriangleMap_->begin();
          iter != nodePtr->externalTriangleMap_->end(); iter++) {
        array<SimplexId, 2> edge = {iter->first.at(1), iter->first.at(2)};
        if(edge[0] > vertexIntervals_[nodePtr->nid - 1]
           && edge[0] <= vertexIntervals_[nodePtr->nid]) {
          (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge) - 1]
            .push_back(iter->second);
        }
      }

      return 0;
    }

    /**
     * Get the triangle edges for all the triangles in a given node.
     */
    int getTriangleEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->triangleEdges_->clear();
      nodePtr->triangleEdges_->resize(
        triangleIntervals_[nodePtr->nid] - triangleIntervals_[nodePtr->nid - 1],
        vector<SimplexId>(3));
      boost::unordered_map<SimplexId, vector<vector<SimplexId>>> edgeNodes;

      if(nodePtr->internalEdgeMap_ == nullptr) {
        nodePtr->internalEdgeMap_
          = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
        buildInternalEdgeMap(
          nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
      }
      if(nodePtr->internalTriangleMap_ == nullptr) {
        nodePtr->internalTriangleMap_
          = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        buildInternalTriangleMap(
          nodePtr->nid, nullptr, nodePtr->internalTriangleMap_, nullptr);
      }

      for(auto iter = nodePtr->internalTriangleMap_->begin();
          iter != nodePtr->internalTriangleMap_->end(); iter++) {
        // since the first vertex of the triangle is in the node ...
        array<SimplexId, 2> edgePair = {iter->first[0], iter->first[1]};
        (*nodePtr->triangleEdges_)[iter->second - 1][0]
          = nodePtr->internalEdgeMap_->at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
        edgePair[1] = iter->first[2];
        (*nodePtr->triangleEdges_)[iter->second - 1][1]
          = nodePtr->internalEdgeMap_->at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
        edgePair[0] = iter->first[1];
        if(edgePair[0] > vertexIntervals_[nodePtr->nid - 1]
           && edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
          (*nodePtr->triangleEdges_)[iter->second - 1][2]
            = nodePtr->internalEdgeMap_->at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
        } else {
          vector<SimplexId> edgeTuple{
            iter->second - 1, edgePair[0], edgePair[1]};
          SimplexId nodeNum = vertexIndices_[edgePair[0]];
          edgeNodes[nodeNum].push_back(edgeTuple);
        }
      }

      ImplicitCluster *exnode;
      for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
        exnode = searchCache(iter->first, nodePtr->nid);
        if(!exnode) {
          boost::unordered_map<array<SimplexId, 2>, SimplexId>
            localInternalEdgeMap;
          buildInternalEdgeMap(
            iter->first, nullptr, &localInternalEdgeMap, nullptr);
          for(vector<SimplexId> edgeTuple : iter->second) {
            array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
            (*nodePtr->triangleEdges_)[edgeTuple[0]][2]
              = localInternalEdgeMap.at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        } else {
          if(exnode->internalEdgeMap_ == nullptr) {
            exnode->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
          }
          for(vector<SimplexId> edgeTuple : iter->second) {
            array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
            (*nodePtr->triangleEdges_)[edgeTuple[0]][2]
              = exnode->internalEdgeMap_->at(edgePair)
                + edgeIntervals_[iter->first - 1];
          }
        }
      }

      return 0;
    }

    /**
     * Get the triangle links for all the triangles in a given node.
     */
    int getTriangleLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->triangleLinks_->clear();
      nodePtr->triangleLinks_->resize(triangleIntervals_[nodePtr->nid]
                                      - triangleIntervals_[nodePtr->nid - 1]);

      if(nodePtr->triangleStars_ == nullptr) {
        if(nodePtr->internalTriangleMap_ == nullptr) {
          nodePtr->internalTriangleMap_
            = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        }
        nodePtr->triangleStars_ = new vector<vector<SimplexId>>();
        buildInternalTriangleMap(nodePtr->nid, nullptr,
                                 nodePtr->internalTriangleMap_,
                                 nodePtr->triangleStars_);
      }

      boost::unordered_map<array<SimplexId, 3>, SimplexId>::const_iterator iter;
      for(iter = nodePtr->internalTriangleMap_->begin();
          iter != nodePtr->internalTriangleMap_->end(); iter++) {
        for(SimplexId i = 0;
            i
            < (SimplexId)(*(nodePtr->triangleStars_))[iter->second - 1].size();
            i++) {
          for(int j = 0; j < 4; j++) {
            SimplexId vertexId = cellArray_->getCellVertex(
              (*(nodePtr->triangleStars_))[iter->second - 1][i], j);
            if((vertexId != iter->first[0]) && (vertexId != iter->first[1])
               && (vertexId != iter->first[2])) {
              (*(nodePtr->triangleLinks_))[iter->second - 1].push_back(
                vertexId);
              break;
            }
          }
        }
      }

      return 0;
    }

    /**
     * Get the vertex edges for all the vertices in a given node.
     */
    int getVertexEdges(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->vertexEdges_->clear();
      nodePtr->vertexEdges_->resize(vertexIntervals_[nodePtr->nid]
                                    - vertexIntervals_[nodePtr->nid - 1]);

      if(nodePtr->internalEdgeMap_ == nullptr) {
        nodePtr->internalEdgeMap_
          = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
        buildInternalEdgeMap(
          nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
      }
      if(nodePtr->externalEdgeMap_ == nullptr) {
        buildExternalEdgeMap(nodePtr);
      }

      boost::unordered_map<array<SimplexId, 2>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalEdgeMap_->begin();
          iter != nodePtr->internalEdgeMap_->end(); iter++) {
        (*(nodePtr->vertexEdges_))[iter->first[0]
                                   - vertexIntervals_[nodePtr->nid - 1] - 1]
          .push_back(edgeIntervals_[nodePtr->nid - 1] + iter->second);
        if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
          (*(nodePtr->vertexEdges_))[iter->first[1]
                                     - vertexIntervals_[nodePtr->nid - 1] - 1]
            .push_back(edgeIntervals_[nodePtr->nid - 1] + iter->second);
        }
      }

      for(iter = nodePtr->externalEdgeMap_->begin();
          iter != nodePtr->externalEdgeMap_->end(); iter++) {
        (*(nodePtr->vertexEdges_))[iter->first[1]
                                   - vertexIntervals_[nodePtr->nid - 1] - 1]
          .push_back(iter->second);
      }

      return 0;
    }

    /**
     * Get the vertex links for all the vertices in a given node.
     */
    int getVertexLinks(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->vertexLinks_->clear();
      nodePtr->vertexLinks_->resize(vertexIntervals_[nodePtr->nid]
                                    - vertexIntervals_[nodePtr->nid - 1]);
      boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;

      if(getDimensionality() == 2) {

        if(nodePtr->internalEdgeMap_ == nullptr) {
          nodePtr->internalEdgeMap_
            = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
          buildInternalEdgeMap(
            nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
        }

        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {
          // the first vertex of the cell must be in the cluster
          array<SimplexId, 2> edgePair
            = {(SimplexId)cellArray_->getCellVertex(cid, 1),
               (SimplexId)cellArray_->getCellVertex(cid, 2)};
          SimplexId nodeId = vertexIndices_[cellArray_->getCellVertex(cid, 1)];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            ImplicitCluster *newNode = new ImplicitCluster(nodeId);
            newNode->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              nodeId, nullptr, newNode->internalEdgeMap_, nullptr);
            nodeMaps[nodeId] = newNode;
          }
          (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 0)
                                     - vertexIntervals_[nodePtr->nid - 1] - 1]
            .push_back(nodeMaps[nodeId]->internalEdgeMap_->at(edgePair)
                       + edgeIntervals_[nodeId - 1]);
          if(cellArray_->getCellVertex(cid, 1)
             <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = cellArray_->getCellVertex(cid, 0);
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 1)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodePtr->internalEdgeMap_->at(edgePair)
                         + edgeIntervals_[nodePtr->nid - 1]);
            if(cellArray_->getCellVertex(cid, 2)
               <= vertexIntervals_[nodePtr->nid]) {
              edgePair[1] = cellArray_->getCellVertex(cid, 1);
              (*(nodePtr
                   ->vertexLinks_))[cellArray_->getCellVertex(cid, 2)
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
                .push_back(nodePtr->internalEdgeMap_->at(edgePair)
                           + edgeIntervals_[nodePtr->nid - 1]);
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          array<SimplexId, 2> edgePair
            = {(SimplexId)cellArray_->getCellVertex(cid, 0),
               (SimplexId)cellArray_->getCellVertex(cid, 2)};
          SimplexId nodeId = vertexIndices_[edgePair[0]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            ImplicitCluster *newNode = new ImplicitCluster(nodeId);
            newNode->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              nodeId, nullptr, newNode->internalEdgeMap_, nullptr);
            nodeMaps[nodeId] = newNode;
          }
          if(cellArray_->getCellVertex(cid, 1)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, 1)
                  <= vertexIntervals_[nodePtr->nid]) {
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 1)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodeMaps[nodeId]->internalEdgeMap_->at(edgePair)
                         + edgeIntervals_[nodeId - 1]);
          }
          if(cellArray_->getCellVertex(cid, 2)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, 2)
                  <= vertexIntervals_[nodePtr->nid]) {
            edgePair[1] = cellArray_->getCellVertex(cid, 1);
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 2)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodeMaps[nodeId]->internalEdgeMap_->at(edgePair)
                         + edgeIntervals_[nodeId - 1]);
          }
        }
      } else if(getDimensionality() == 3) {

        if(nodePtr->internalTriangleMap_ == nullptr) {
          nodePtr->internalTriangleMap_
            = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
          buildInternalTriangleMap(
            nodePtr->nid, nullptr, nodePtr->internalTriangleMap_, nullptr);
        }

        for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
            cid <= cellIntervals_[nodePtr->nid]; cid++) {

          // v1: (v2, v3, v4)
          array<SimplexId, 3> triangleVec;
          triangleVec[0] = cellArray_->getCellVertex(cid, 1);
          triangleVec[1] = cellArray_->getCellVertex(cid, 2);
          triangleVec[2] = cellArray_->getCellVertex(cid, 3);
          SimplexId nodeId = vertexIndices_[cellArray_->getCellVertex(cid, 1)];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            ImplicitCluster *newNode = new ImplicitCluster(nodeId);
            newNode->internalTriangleMap_
              = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            buildInternalTriangleMap(
              nodeId, nullptr, newNode->internalTriangleMap_, nullptr);
            nodeMaps[nodeId] = newNode;
          }
          (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 0)
                                     - vertexIntervals_[nodePtr->nid - 1] - 1]
            .push_back(nodeMaps[nodeId]->internalTriangleMap_->at(triangleVec)
                       + triangleIntervals_[nodeId - 1]);
          // v2: (v1, v3, v4)
          if(cellArray_->getCellVertex(cid, 1)
             <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[0] = cellArray_->getCellVertex(cid, 0);
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 1)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodePtr->internalTriangleMap_->at(triangleVec)
                         + triangleIntervals_[nodePtr->nid - 1]);
            // v3: (v1, v2, v4)
            if(cellArray_->getCellVertex(cid, 2)
               <= vertexIntervals_[nodePtr->nid]) {
              triangleVec[1] = cellArray_->getCellVertex(cid, 1);
              (*(nodePtr
                   ->vertexLinks_))[cellArray_->getCellVertex(cid, 2)
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
                .push_back(nodePtr->internalTriangleMap_->at(triangleVec)
                           + triangleIntervals_[nodePtr->nid - 1]);
            }
            // v4: (v1, v2, v3)
            if(cellArray_->getCellVertex(cid, 3)
               <= vertexIntervals_[nodePtr->nid]) {
              triangleVec[2] = cellArray_->getCellVertex(cid, 2);
              (*(nodePtr
                   ->vertexLinks_))[cellArray_->getCellVertex(cid, 3)
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
                .push_back(nodePtr->internalTriangleMap_->at(triangleVec)
                           + triangleIntervals_[nodePtr->nid - 1]);
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]) {
          // start from v2
          array<SimplexId, 3> triangleVec;
          triangleVec[0] = cellArray_->getCellVertex(cid, 0);
          triangleVec[1] = cellArray_->getCellVertex(cid, 2);
          triangleVec[2] = cellArray_->getCellVertex(cid, 3);
          SimplexId nodeId = vertexIndices_[triangleVec[0]];
          if(nodeMaps.find(nodeId) == nodeMaps.end()) {
            ImplicitCluster *newNode = new ImplicitCluster(nodeId);
            newNode->internalTriangleMap_
              = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            buildInternalTriangleMap(
              nodeId, nullptr, newNode->internalTriangleMap_, nullptr);
            nodeMaps[nodeId] = newNode;
          }
          if(cellArray_->getCellVertex(cid, 1)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, 1)
                  <= vertexIntervals_[nodePtr->nid]) {
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 1)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodeMaps[nodeId]->internalTriangleMap_->at(triangleVec)
                         + triangleIntervals_[nodeId - 1]);
          }
          if(cellArray_->getCellVertex(cid, 2)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, 2)
                  <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[1] = cellArray_->getCellVertex(cid, 1);
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 2)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodeMaps[nodeId]->internalTriangleMap_->at(triangleVec)
                         + triangleIntervals_[nodeId - 1]);
          }
          if(cellArray_->getCellVertex(cid, 3)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, 3)
                  <= vertexIntervals_[nodePtr->nid]) {
            triangleVec[1] = cellArray_->getCellVertex(cid, 1);
            triangleVec[2] = cellArray_->getCellVertex(cid, 2);
            (*(nodePtr->vertexLinks_))[cellArray_->getCellVertex(cid, 3)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(nodeMaps[nodeId]->internalTriangleMap_->at(triangleVec)
                         + triangleIntervals_[nodeId - 1]);
          }
        }
      }

      // release the memory
      for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
        delete iter->second;
      }

      return 0;
    }

    /**
     * Get the vertex neighbors for all the vertices in a given node.
     */
    int getVertexNeighbors(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->vertexNeighbors_->clear();
      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
      SimplexId localVertexNum
        = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
      nodePtr->vertexNeighbors_->resize(localVertexNum);

      SimplexId v1, v2;
      vector<boost::unordered_set<SimplexId>> vertexNeighborSet(localVertexNum);

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

      for(SimplexId i = 0; i < localVertexNum; i++) {
        (*(nodePtr->vertexNeighbors_))[i] = vector<SimplexId>(
          vertexNeighborSet[i].begin(), vertexNeighborSet[i].end());
      }

      return 0;
    }

    /**
     * Get the vertex stars for all the vertices in a given node.
     * The function is similar as getVertexCells().
     */
    int getVertexStars(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->vertexStars_->clear();
      nodePtr->vertexStars_->resize(vertexIntervals_[nodePtr->nid]
                                    - vertexIntervals_[nodePtr->nid - 1]);
      SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);

      // loop through the internal cell list
      for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
          cid <= cellIntervals_[nodePtr->nid]; cid++) {
        for(SimplexId j = 0; j < verticesPerCell; j++) {
          // see if it is in the current node
          if(cellArray_->getCellVertex(cid, j)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, j)
                  <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr->vertexStars_))[cellArray_->getCellVertex(cid, j)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(cid);
        }
      }

      // and also external cell list
      for(SimplexId cid : externalCells_[nodePtr->nid]) {
        for(SimplexId j = 0; j < verticesPerCell; j++) {
          // see if it is in the current node
          if(cellArray_->getCellVertex(cid, j)
               > vertexIntervals_[nodePtr->nid - 1]
             && cellArray_->getCellVertex(cid, j)
                  <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr->vertexStars_))[cellArray_->getCellVertex(cid, j)
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(cid);
        }
      }

      return 0;
    }

    /**
     * Get the vertex triangles for all the vertices in a given node.
     */
    int getVertexTriangles(ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
        return -1;
#endif

      nodePtr->vertexTriangles_->clear();
      nodePtr->vertexTriangles_->resize(vertexIntervals_[nodePtr->nid]
                                        - vertexIntervals_[nodePtr->nid - 1]);

      if(nodePtr->internalTriangleMap_ == nullptr) {
        nodePtr->internalTriangleMap_
          = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
        buildInternalTriangleMap(
          nodePtr->nid, nullptr, nodePtr->internalTriangleMap_, nullptr);
      }

      if(nodePtr->externalTriangleMap_ == nullptr) {
        buildExternalTriangleMap(nodePtr);
      }

      boost::unordered_map<array<SimplexId, 3>, SimplexId>::iterator iter;
      for(iter = nodePtr->internalTriangleMap_->begin();
          iter != nodePtr->internalTriangleMap_->end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[j] <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr
                 ->vertexTriangles_))[iter->first[j]
                                      - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(iter->second + triangleIntervals_[nodePtr->nid - 1]);
        }
      }

      for(iter = nodePtr->externalTriangleMap_->begin();
          iter != nodePtr->externalTriangleMap_->end(); iter++) {
        for(SimplexId j = 0; j < 3; j++) {
          if(iter->first.at(j) > vertexIntervals_[nodePtr->nid - 1]
             && iter->first.at(j) <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr
                 ->vertexTriangles_))[iter->first.at(j)
                                      - vertexIntervals_[nodePtr->nid - 1] - 1]
              .push_back(iter->second);
        }
      }

      return 0;
    }

    /**
     * Get the boundary cells in a given node.
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
        if(nodePtr->boundaryEdges_ == nullptr) {
          nodePtr->boundaryEdges_ = new vector<bool>(localEdgeNum, false);
          if(nodePtr->edgeStars_ == nullptr) {
            if(nodePtr->internalEdgeMap_ == nullptr) {
              nodePtr->internalEdgeMap_
                = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            }
            nodePtr->edgeStars_ = new vector<vector<SimplexId>>();
            buildInternalEdgeMap(nodePtr->nid, nullptr,
                                 nodePtr->internalEdgeMap_,
                                 nodePtr->edgeStars_);
          }
          for(SimplexId i = 0; i < localEdgeNum; i++) {
            if((*(nodePtr->edgeStars_))[i].size() == 1) {
              (*(nodePtr->boundaryEdges_))[i] = true;
            }
          }
        }
        // if boundary vertices are requested
        if(dim == 0 && nodePtr->boundaryVertices_ == nullptr) {
          nodePtr->boundaryVertices_ = new vector<bool>(
            vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalEdgeMap_ == nullptr) {
            nodePtr->externalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildExternalEdgeMap(nodePtr);
          }
          // internal edges
          for(auto iter = nodePtr->internalEdgeMap_->begin();
              iter != nodePtr->internalEdgeMap_->end(); iter++) {
            if((*(nodePtr->boundaryEdges_))[iter->second - 1]) {
              (*(nodePtr->boundaryVertices_))
                [iter->first[0] - vertexIntervals_[nodePtr->nid - 1] - 1]
                = true;
              if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                (*(nodePtr->boundaryVertices_))
                  [iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
            }
          }
          // external edges
          boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;
          for(auto iter = nodePtr->externalEdgeMap_->begin();
              iter != nodePtr->externalEdgeMap_->end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = new ImplicitCluster(nodeId);
              getBoundaryCells(nodeMaps[nodeId]);
            }
            if((*(nodeMaps[nodeId]
                    ->boundaryEdges_))[iter->second - edgeIntervals_[nodeId - 1]
                                       - 1]) {
              (*(nodePtr->boundaryVertices_))
                [iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1]
                = true;
            }
          }

          // release the memory
          for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
            delete iter->second;
          }
        }
      } else if(getDimensionality() == 3) {
        // get the boundary triangles first
        SimplexId localTriangleNum = triangleIntervals_[nodePtr->nid]
                                     - triangleIntervals_[nodePtr->nid - 1];
        if(nodePtr->boundaryTriangles_ == nullptr) {
          nodePtr->boundaryTriangles_
            = new vector<bool>(localTriangleNum, false);
          if(nodePtr->triangleStars_ == nullptr) {
            if(nodePtr->internalTriangleMap_ == nullptr) {
              nodePtr->internalTriangleMap_
                = new boost::unordered_map<array<SimplexId, 3>, SimplexId>();
            }
            nodePtr->triangleStars_ = new vector<vector<SimplexId>>();
            buildInternalTriangleMap(nodePtr->nid, nullptr,
                                     nodePtr->internalTriangleMap_,
                                     nodePtr->triangleStars_);
          }
          for(SimplexId i = 0; i < localTriangleNum; i++) {
            if((*(nodePtr->triangleStars_))[i].size() == 1) {
              (*(nodePtr->boundaryTriangles_))[i] = true;
            }
          }
        }
        // if the boundary edges are requested
        if(dim == 1 && nodePtr->boundaryEdges_ == nullptr) {
          nodePtr->boundaryEdges_ = new vector<bool>(
            edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalTriangleMap_ == nullptr) {
            buildExternalTriangleMap(nodePtr);
          }
          if(nodePtr->internalEdgeMap_ == nullptr) {
            nodePtr->internalEdgeMap_
              = new boost::unordered_map<array<SimplexId, 2>, SimplexId>();
            buildInternalEdgeMap(
              nodePtr->nid, nullptr, nodePtr->internalEdgeMap_, nullptr);
          }
          // internal triangles
          for(auto iter = nodePtr->internalTriangleMap_->begin();
              iter != nodePtr->internalTriangleMap_->end(); iter++) {
            if((*(nodePtr->boundaryTriangles_))[iter->second - 1]) {
              array<SimplexId, 2> edgePair = {iter->first[0], iter->first[1]};
              (*(nodePtr
                   ->boundaryEdges_))[nodePtr->internalEdgeMap_->at(edgePair)
                                      - 1]
                = true;
              edgePair[1] = iter->first[2];
              (*(nodePtr
                   ->boundaryEdges_))[nodePtr->internalEdgeMap_->at(edgePair)
                                      - 1]
                = true;
              if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                edgePair[0] = iter->first[1];
                (*(nodePtr
                     ->boundaryEdges_))[nodePtr->internalEdgeMap_->at(edgePair)
                                        - 1]
                  = true;
              }
            }
          }
          // external triangles
          boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;
          for(auto iter = nodePtr->externalTriangleMap_->begin();
              iter != nodePtr->externalTriangleMap_->end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = new ImplicitCluster(nodeId);
              getBoundaryCells(nodeMaps[nodeId]);
            }
            if((*(nodeMaps[nodeId]->boundaryTriangles_))
                 [iter->second - triangleIntervals_[nodeId - 1] - 1]) {
              if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                array<SimplexId, 2> edgePair = {iter->first[1], iter->first[2]};
                (*(nodePtr
                     ->boundaryEdges_))[nodePtr->internalEdgeMap_->at(edgePair)
                                        - 1]
                  = true;
              }
            }
          }

          // release the memory
          for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
            delete iter->second;
          }
        }

        // if the boundary vertices are requested
        else if(dim == 0 && nodePtr->boundaryVertices_ == nullptr) {
          nodePtr->boundaryVertices_ = new vector<bool>(
            vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
            false);
          if(nodePtr->externalTriangleMap_ == nullptr) {
            buildExternalTriangleMap(nodePtr);
          }
          // internal triangles
          for(auto iter = nodePtr->internalTriangleMap_->begin();
              iter != nodePtr->internalTriangleMap_->end(); iter++) {
            if((*(nodePtr->boundaryTriangles_))[iter->second - 1]) {
              for(int j = 0; j < 3; j++) {
                SimplexId vid = iter->first[j];
                if(vid <= vertexIntervals_[nodePtr->nid]) {
                  (*(nodePtr->boundaryVertices_))
                    [vid - vertexIntervals_[nodePtr->nid - 1] - 1]
                    = true;
                }
              }
            }
          }
          // external triangles
          boost::unordered_map<SimplexId, ImplicitCluster *> nodeMaps;
          for(auto iter = nodePtr->externalTriangleMap_->begin();
              iter != nodePtr->externalTriangleMap_->end(); iter++) {
            SimplexId nodeId = vertexIndices_[iter->first[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = new ImplicitCluster(nodeId);
              getBoundaryCells(nodeMaps[nodeId]);
            }
            if((*(nodeMaps[nodeId]->boundaryTriangles_))
                 [iter->second - triangleIntervals_[nodeId - 1] - 1]) {
              if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
                (*(nodePtr->boundaryVertices_))
                  [iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
              if(iter->first[2] > vertexIntervals_[nodePtr->nid - 1]
                 && iter->first[2] <= vertexIntervals_[nodePtr->nid]) {
                (*(nodePtr->boundaryVertices_))
                  [iter->first[2] - vertexIntervals_[nodePtr->nid - 1] - 1]
                  = true;
              }
            }
          }

          // release the memory
          for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++) {
            delete iter->second;
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
    vector<SimplexId> vertexIntervals_;
    vector<SimplexId> edgeIntervals_;
    vector<SimplexId> triangleIntervals_;
    vector<SimplexId> cellIntervals_;
    std::shared_ptr<CellArray> cellArray_;
    vector<vector<SimplexId>> externalCells_;

    // Cache system
    size_t cacheSize_;
    mutable vector<list<ImplicitCluster *>> caches_;
    mutable vector<
      boost::unordered_map<SimplexId, list<ImplicitCluster *>::iterator>>
      cacheMaps_;

    friend class TestTopoCluster;
  };
} // namespace ttk