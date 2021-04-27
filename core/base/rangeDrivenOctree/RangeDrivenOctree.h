/// \ingroup base
/// \class ttk::RangeDrivenOctree
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2015.
///
/// \brief TTK optional package for octree based range queries in bivariate
/// volumetric data.
///
/// This class accelerates range-driven queries in bivariate volumetric data.
/// This class is typically used to accelerate fiber surface computation.
///
/// \b Related \b publication \n
/// "Fast and Exact Fiber Surface Extraction for Tetrahedral Meshes" \n
/// Pavol Klacansky, Julien Tierny, Hamish Carr, Zhao Geng \n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa FiberSurface.h %for a usage example.

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  class RangeDrivenOctree : virtual public Debug {
  public:
    RangeDrivenOctree();

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int build(const triangulationType *const triangulation);

    inline bool empty() const {
      return nodeList_.empty();
    }

    void flush();

    int getTet2NodeMap(std::vector<SimplexId> &map,
                       const bool &forSegmentation = false) const;

    int rangeSegmentQuery(const std::pair<double, double> &p0,
                          const std::pair<double, double> &p1,
                          std::vector<SimplexId> &cellList) const;

    inline void setCellList(const SimplexId *cellList) {
      cellList_ = cellList;
    }

    inline void setCellNumber(const SimplexId &cellNumber) {
      cellNumber_ = cellNumber;
    }

    inline void setLeafMinimumCellNumber(const SimplexId &number) {
      leafMinimumCellNumber_ = number;
    }

    inline void setLeafMinimumDomainVolumeRatio(const float &ratio) {
      leafMinimumRangeAreaRatio_ = ratio;
    }

    inline void setLeafMinimumRangeAreaRatio(const float &ratio) {
      leafMinimumRangeAreaRatio_ = ratio;
    }

    inline void setPointList(const float *pointList) {
      pointList_ = pointList;
    }

    inline void setRange(const void *u, const void *v) {
      u_ = u;
      v_ = v;
    }

    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    int stats(std::ostream &stream);

    int statNode(const SimplexId &nodeId, std::ostream &stream);

  protected:
    struct OctreeNode {
      std::pair<std::pair<double, double>, std::pair<double, double>>
        rangeBox_{};
      std::vector<SimplexId> cellList_{};
      std::vector<SimplexId> childList_{};
      std::vector<std::pair<float, float>> domainBox_{};
    };

    template <class dataTypeU, class dataTypeV>
    int buildNode(const std::vector<SimplexId> &cellList,
                  const std::vector<std::pair<float, float>> &domainBox,
                  const std::pair<std::pair<double, double>,
                                  std::pair<double, double>> &rangeBox,
                  SimplexId &nodeId);

    int rangeSegmentQuery(const std::pair<double, double> &p0,
                          const std::pair<double, double> &p1,
                          const SimplexId &nodeId,
                          std::vector<SimplexId> &cellList) const;

    bool segmentIntersection(const std::pair<double, double> &p0,
                             const std::pair<double, double> &p1,
                             const std::pair<double, double> &q0,
                             const std::pair<double, double> &q1) const;

    const void *u_{};
    const void *v_{};
    const float *pointList_{};
    const SimplexId *cellList_{};
    float domainVolume_{}, leafMinimumDomainVolumeRatio_{0.01F},
      leafMinimumRangeAreaRatio_{0.01F}, rangeArea_{};
    SimplexId cellNumber_{}, vertexNumber_{}, leafMinimumCellNumber_{6},
      rootId_{};
    mutable SimplexId queryResultNumber_{};
    std::vector<OctreeNode> nodeList_{};
    std::vector<std::vector<std::pair<float, float>>> cellDomainBox_{};
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
      cellRangeBox_{};
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <RangeDrivenOctree.cpp>

template <class dataTypeU, class dataTypeV, typename triangulationType>
int ttk::RangeDrivenOctree::build(
  const triangulationType *const triangulation) {

  Timer t;

  dataTypeU *u = (dataTypeU *)u_;
  dataTypeV *v = (dataTypeV *)v_;

  if(triangulation) {
    cellNumber_ = triangulation->getNumberOfCells();
  }

  if(triangulation) {
    vertexNumber_ = triangulation->getNumberOfVertices();
  }

  cellDomainBox_.resize(cellNumber_, std::vector<std::pair<float, float>>(3));

  cellRangeBox_.resize(cellNumber_);

  // WARNING: assuming tets only here
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber_; i++) {

    for(int j = 0; j < 3; j++) {
      cellDomainBox_[i][j].first = FLT_MAX;
      cellDomainBox_[i][j].second = -FLT_MAX;
    }

    const SimplexId *cell = NULL;
    if(!triangulation) {
      cell = &(cellList_[5 * i + 1]);
    }

    for(int j = 0; j < 4; j++) {

      // update the domain bounding box for that tet
      SimplexId vertexId = 0;

      if(triangulation) {
        triangulation->getCellVertex(i, j, vertexId);
      } else if(cell != nullptr) {
        vertexId = cell[j];
      }

      std::array<float, 3> p{};
      if(triangulation) {
        triangulation->getVertexPoint(vertexId, p[0], p[1], p[2]);
      } else {
        p[0] = pointList_[3 * vertexId];
        p[1] = pointList_[3 * vertexId + 1];
        p[2] = pointList_[3 * vertexId + 2];
      }

      // update x-range
      if(p[0] < cellDomainBox_[i][0].first)
        cellDomainBox_[i][0].first = p[0];
      if(p[0] > cellDomainBox_[i][0].second)
        cellDomainBox_[i][0].second = p[0];

      // update y-range
      if(p[1] < cellDomainBox_[i][1].first)
        cellDomainBox_[i][1].first = p[1];
      if(p[1] > cellDomainBox_[i][1].second)
        cellDomainBox_[i][1].second = p[1];

      // update z-range
      if(p[2] < cellDomainBox_[i][2].first)
        cellDomainBox_[i][2].first = p[2];
      if(p[2] > cellDomainBox_[i][2].second)
        cellDomainBox_[i][2].second = p[2];

      // update the range bounding box
      if(!j) {
        cellRangeBox_[i].first.first = u[vertexId];
        cellRangeBox_[i].first.second = u[vertexId];

        cellRangeBox_[i].second.first = v[vertexId];
        cellRangeBox_[i].second.second = v[vertexId];
      } else {

        // update u
        if(u[vertexId] < cellRangeBox_[i].first.first)
          cellRangeBox_[i].first.first = u[vertexId];
        if(u[vertexId] > cellRangeBox_[i].first.second)
          cellRangeBox_[i].first.second = u[vertexId];

        // update v
        if(v[vertexId] < cellRangeBox_[i].second.first)
          cellRangeBox_[i].second.first = v[vertexId];
        if(v[vertexId] > cellRangeBox_[i].second.second)
          cellRangeBox_[i].second.second = v[vertexId];
      }
    }
  }

  std::vector<SimplexId> domain(cellNumber_);
  for(SimplexId i = 0; i < cellNumber_; i++)
    domain[i] = i;

  // get global bBoxes
  std::vector<std::pair<float, float>> domainBox(3);
  std::pair<std::pair<double, double>, std::pair<double, double>> rangeBox;

  for(SimplexId i = 0; i < vertexNumber_; i++) {

    // domain one
    std::array<float, 3> p{};
    if(triangulation) {
      triangulation->getVertexPoint(i, p[0], p[1], p[2]);
    } else {
      p[0] = pointList_[3 * i];
      p[1] = pointList_[3 * i + 1];
      p[2] = pointList_[3 * i + 2];
    }

    for(int j = 0; j < 3; j++) {
      if(!i) {
        domainBox[j].first = domainBox[j].second = p[j];
      } else {
        if(p[j] < domainBox[j].first)
          domainBox[j].first = p[j];
        if(p[j] > domainBox[j].second)
          domainBox[j].second = p[j];
      }
    }

    if(!i) {
      rangeBox.first.first = rangeBox.first.second = u[i];
      rangeBox.second.first = rangeBox.second.second = v[i];
    } else {
      if(u[i] < rangeBox.first.first)
        rangeBox.first.first = u[i];
      if(u[i] > rangeBox.first.second)
        rangeBox.first.second = u[i];

      if(v[i] < rangeBox.second.first)
        rangeBox.second.first = v[i];
      if(v[i] > rangeBox.second.second)
        rangeBox.second.second = v[i];
    }
  }

  rangeArea_ = (rangeBox.first.second - rangeBox.first.first)
               * (rangeBox.second.second - rangeBox.second.first);
  domainVolume_ = (domainBox[0].second - domainBox[0].first)
                  * (domainBox[1].second - domainBox[1].first)
                  * (domainBox[2].second - domainBox[2].first);

  // special case for tets obtained from regular grid subdivision (assuming 6)
  if(leafMinimumCellNumber_ < 6)
    leafMinimumCellNumber_ = 6;

  leafMinimumDomainVolumeRatio_ = (1.0 / ((float)cellNumber_)) / 2.0;

  this->printMsg(
    "Range area ratio: " + std::to_string(leafMinimumRangeAreaRatio_),
    debug::Priority::DETAIL);

  buildNode<dataTypeU, dataTypeV>(domain, domainBox, rangeBox, rootId_);

  this->printMsg("Octree built", 1.0, t.getElapsedTime(), this->threadNumber_);

  // debug
  //   stats(std::cout);
  // end of debug

  return 0;
}

template <class dataTypeU, class dataTypeV>
int ttk::RangeDrivenOctree::buildNode(
  const std::vector<SimplexId> &cellList,
  const std::vector<std::pair<float, float>> &domainBox,
  const std::pair<std::pair<double, double>, std::pair<double, double>>
    &rangeBox,
  SimplexId &nodeId) {

  nodeId = nodeList_.size();

  nodeList_.resize(nodeList_.size() + 1);

  nodeList_.back().rangeBox_ = rangeBox;
  nodeList_.back().domainBox_ = domainBox;

  float rangeArea = (rangeBox.first.second - rangeBox.first.first)
                    * (rangeBox.second.second - rangeBox.second.first);

  float domainVolume = (domainBox[0].second - domainBox[0].first)
                       * (domainBox[1].second - domainBox[1].first)
                       * (domainBox[2].second - domainBox[2].first);

  if(((SimplexId)cellList.size() > leafMinimumCellNumber_)
     && (rangeArea > leafMinimumRangeAreaRatio_ * rangeArea_)
     && (domainVolume > leafMinimumDomainVolumeRatio_ * domainVolume_)) {

    nodeList_.back().childList_.resize(8);

    std::vector<std::vector<std::pair<float, float>>> childDomainBox(8);
    for(SimplexId i = 0; i < (SimplexId)childDomainBox.size(); i++) {
      childDomainBox[i].resize(3);
    }
    std::vector<std::pair<std::pair<dataTypeU, dataTypeU>,
                          std::pair<dataTypeV, dataTypeV>>>
      childRangeBox(8);
    std::vector<std::vector<SimplexId>> childCellList(8);

    float midX
      = domainBox[0].first + (domainBox[0].second - domainBox[0].first) / 2.0;
    float midY
      = domainBox[1].first + (domainBox[1].second - domainBox[1].first) / 2.0;
    float midZ
      = domainBox[2].first + (domainBox[2].second - domainBox[2].first) / 2.0;

    // 0 - - -
    childDomainBox[0][0].first = domainBox[0].first;
    childDomainBox[0][0].second = midX;
    childDomainBox[0][1].first = domainBox[1].first;
    childDomainBox[0][1].second = midY;
    childDomainBox[0][2].first = domainBox[2].first;
    childDomainBox[0][2].second = midZ;

    // 1 - - +
    childDomainBox[1][0].first = domainBox[0].first;
    childDomainBox[1][0].second = midX;
    childDomainBox[1][1].first = domainBox[1].first;
    childDomainBox[1][1].second = midY;
    childDomainBox[1][2].first = midZ;
    childDomainBox[1][2].second = domainBox[2].second;

    // 2 - + -
    childDomainBox[2][0].first = domainBox[0].first;
    childDomainBox[2][0].second = midX;
    childDomainBox[2][1].first = midY;
    childDomainBox[2][1].second = domainBox[1].second;
    childDomainBox[2][2].first = domainBox[2].first;
    childDomainBox[2][2].second = midZ;

    // 3 - + +
    childDomainBox[3][0].first = domainBox[0].first;
    childDomainBox[3][0].second = midX;
    childDomainBox[3][1].first = midY;
    childDomainBox[3][1].second = domainBox[1].second;
    childDomainBox[3][2].first = midZ;
    childDomainBox[3][2].second = domainBox[2].second;

    // 4 + - -
    childDomainBox[4][0].first = midX;
    childDomainBox[4][0].second = domainBox[0].second;
    childDomainBox[4][1].first = domainBox[1].first;
    childDomainBox[4][1].second = midY;
    childDomainBox[4][2].first = domainBox[2].first;
    childDomainBox[4][2].second = midZ;

    // 5 + - +
    childDomainBox[5][0].first = midX;
    childDomainBox[5][0].second = domainBox[0].second;
    childDomainBox[5][1].first = domainBox[1].first;
    childDomainBox[5][1].second = midY;
    childDomainBox[5][2].first = midZ;
    childDomainBox[5][2].second = domainBox[2].second;

    // 6 + + -
    childDomainBox[6][0].first = midX;
    childDomainBox[6][0].second = domainBox[0].second;
    childDomainBox[6][1].first = midY;
    childDomainBox[6][1].second = domainBox[1].second;
    childDomainBox[6][2].first = domainBox[2].first;
    childDomainBox[6][2].second = midZ;

    // 7 + + +
    childDomainBox[7][0].first = midX;
    childDomainBox[7][0].second = domainBox[0].second;
    childDomainBox[7][1].first = midY;
    childDomainBox[7][1].second = domainBox[1].second;
    childDomainBox[7][2].first = midZ;
    childDomainBox[7][2].second = domainBox[2].second;

    for(SimplexId i = 0; i < (SimplexId)cellList.size(); i++) {

      SimplexId childId = 0;

      for(int j = 0; j < 8; j++) {
        if((cellDomainBox_[cellList[i]][0].first >= childDomainBox[j][0].first)
           && (cellDomainBox_[cellList[i]][0].first
               < childDomainBox[j][0].second)
           && (cellDomainBox_[cellList[i]][1].first
               >= childDomainBox[j][1].first)
           && (cellDomainBox_[cellList[i]][1].first
               < childDomainBox[j][1].second)
           && (cellDomainBox_[cellList[i]][2].first
               >= childDomainBox[j][2].first)
           && (cellDomainBox_[cellList[i]][2].first
               < childDomainBox[j][2].second)) {

          childId = j;
          break;
        }
      }

      // update child's range box
      if(childCellList[childId].empty()) {
        childRangeBox[childId].first.first
          = cellRangeBox_[cellList[i]].first.first;
        childRangeBox[childId].first.second
          = cellRangeBox_[cellList[i]].first.second;

        childRangeBox[childId].second.first
          = cellRangeBox_[cellList[i]].second.first;
        childRangeBox[childId].second.second
          = cellRangeBox_[cellList[i]].second.second;
      } else {
        if(cellRangeBox_[cellList[i]].first.first
           < childRangeBox[childId].first.first) {
          childRangeBox[childId].first.first
            = cellRangeBox_[cellList[i]].first.first;
        }
        if(cellRangeBox_[cellList[i]].first.second
           > childRangeBox[childId].first.second) {
          childRangeBox[childId].first.second
            = cellRangeBox_[cellList[i]].first.second;
        }

        if(cellRangeBox_[cellList[i]].second.first
           < childRangeBox[childId].second.first) {
          childRangeBox[childId].second.first
            = cellRangeBox_[cellList[i]].second.first;
        }
        if(cellRangeBox_[cellList[i]].second.second
           > childRangeBox[childId].second.second) {
          childRangeBox[childId].second.second
            = cellRangeBox_[cellList[i]].second.second;
        }
      }

      childCellList[childId].push_back(cellList[i]);
    }

    for(int i = 0; i < 8; i++) {
      buildNode<dataTypeU, dataTypeV>(childCellList[i], childDomainBox[i],
                                      childRangeBox[i],
                                      nodeList_[nodeId].childList_[i]);
    }
  } else {
    // leaf
    nodeList_[nodeId].cellList_ = cellList;
  }

  return 0;
}
