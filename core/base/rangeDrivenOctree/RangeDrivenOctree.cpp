// NOTE
// this class has been implemented recursively. this is bad programming.
// a non-recursive version will be implemented asap.

#include <RangeDrivenOctree.h>

using namespace ttk;

RangeDrivenOctree::RangeDrivenOctree() {
  this->setDebugMsgPrefix("RangeDrivenOctree");
}

void RangeDrivenOctree::flush() {
  nodeList_.clear();
  cellDomainBox_.clear();
  cellRangeBox_.clear();
}

int RangeDrivenOctree::getTet2NodeMap(std::vector<SimplexId> &map,
                                      const bool &forSegmentation) const {

  std::vector<SimplexId> randomMap;
  if(forSegmentation) {
    randomMap.resize(nodeList_.size());
    for(size_t i = 0; i < randomMap.size(); i++) {
      randomMap[i] = rand() % (randomMap.size());
    }
  }

  map.resize(cellNumber_);
  for(size_t i = 0; i < nodeList_.size(); i++) {
    for(size_t j = 0; j < nodeList_[i].cellList_.size(); j++) {
      if(forSegmentation) {
        map[nodeList_[i].cellList_[j]] = randomMap[i];
      } else {
        map[nodeList_[i].cellList_[j]] = i;
      }
    }
  }

  return 0;
}

int RangeDrivenOctree::rangeSegmentQuery(
  const std::pair<double, double> &p0,
  const std::pair<double, double> &p1,
  std::vector<SimplexId> &cellList) const {

  Timer t;

  queryResultNumber_ = 0;
  cellList.clear();

  SimplexId ret = rangeSegmentQuery(p0, p1, rootId_, cellList);

  this->printMsg("Query done", 1.0, t.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);
  this->printMsg(
    std::vector<std::vector<std::string>>{
      {"#Non empty leaves", std::to_string(cellList.size())}},
    debug::Priority::DETAIL);

  return ret;
}

int RangeDrivenOctree::rangeSegmentQuery(
  const std::pair<double, double> &p0,
  const std::pair<double, double> &p1,
  const SimplexId &nodeId,
  std::vector<SimplexId> &cellList) const {

  // check for intersection for each segment of the range bounding box
  std::pair<double, double> q0, q1;

  // bottom range segment (min, min) (max, min)
  q0.first = nodeList_[nodeId].rangeBox_.first.first;
  q0.second = nodeList_[nodeId].rangeBox_.second.first;

  q1.first = nodeList_[nodeId].rangeBox_.first.second;
  q1.second = q0.second;

  if(segmentIntersection(p0, p1, q0, q1)) {
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }

  // right segment (max, min) (max, max)
  q0.first = nodeList_[nodeId].rangeBox_.first.second;
  q0.second = nodeList_[nodeId].rangeBox_.second.first;

  q1.first = nodeList_[nodeId].rangeBox_.first.second;
  q1.second = nodeList_[nodeId].rangeBox_.second.second;

  if(segmentIntersection(p0, p1, q0, q1)) {
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }

  // top segment (min, max) (max, max)
  q0.first = nodeList_[nodeId].rangeBox_.first.first;
  q0.second = nodeList_[nodeId].rangeBox_.second.second;

  q1.first = nodeList_[nodeId].rangeBox_.first.second;
  q1.second = nodeList_[nodeId].rangeBox_.second.second;

  if(segmentIntersection(p0, p1, q0, q1)) {
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }

  // left segment (min, min) (min, max)
  q0.first = nodeList_[nodeId].rangeBox_.first.first;
  q0.second = nodeList_[nodeId].rangeBox_.second.first;

  q1.first = nodeList_[nodeId].rangeBox_.first.first;
  q1.second = nodeList_[nodeId].rangeBox_.second.second;

  if(segmentIntersection(p0, p1, q0, q1)) {
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }

  // is the segment completely included in the range bounding box?
  if((p0.first >= nodeList_[nodeId].rangeBox_.first.first)
     && (p0.first < nodeList_[nodeId].rangeBox_.first.second)
     && (p0.second >= nodeList_[nodeId].rangeBox_.second.first)
     && (p0.second < nodeList_[nodeId].rangeBox_.second.second)) {

    // p0 is in there
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  if((p1.first >= nodeList_[nodeId].rangeBox_.first.first)
     && (p1.first < nodeList_[nodeId].rangeBox_.first.second)
     && (p1.second >= nodeList_[nodeId].rangeBox_.second.first)
     && (p1.second < nodeList_[nodeId].rangeBox_.second.second)) {

    // p1 is in there
    if(nodeList_[nodeId].childList_.size()) {
      for(size_t i = 0; i < nodeList_[nodeId].childList_.size(); i++) {
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    } else {
      // terminal leaf
      // return our cells
      this->printMsg("Node #" + std::to_string(nodeId) + " returns its "
                       + std::to_string(nodeList_[nodeId].cellList_.size())
                       + " cells(s).",
                     debug::Priority::VERBOSE);

      cellList.insert(cellList.end(), nodeList_[nodeId].cellList_.begin(),
                      nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }

  return 0;
}

int RangeDrivenOctree::statNode(const SimplexId &nodeId, std::ostream &stream) {

  stream << "[RangeDrivenOctree]" << std::endl;
  stream << "[RangeDrivenOctree] Node #" << nodeId << std::endl;
  stream << "[RangeDrivenOctree]   Domain box: ["
         << nodeList_[nodeId].domainBox_[0].first << " "
         << nodeList_[nodeId].domainBox_[0].second << "] ["
         << nodeList_[nodeId].domainBox_[1].first << " "
         << nodeList_[nodeId].domainBox_[1].second << "] ["
         << nodeList_[nodeId].domainBox_[2].first << " "
         << nodeList_[nodeId].domainBox_[2].second << "] "
         << " volume="
         << (nodeList_[nodeId].domainBox_[0].second
             - nodeList_[nodeId].domainBox_[0].first)
              * (nodeList_[nodeId].domainBox_[1].second
                 - nodeList_[nodeId].domainBox_[1].first)
              * (nodeList_[nodeId].domainBox_[2].second
                 - nodeList_[nodeId].domainBox_[2].first)
         << " threshold=" << leafMinimumDomainVolumeRatio_ * domainVolume_
         << std::endl;
  stream << "[RangeDrivenOctree]   Range box: ["
         << nodeList_[nodeId].rangeBox_.first.first << " "
         << nodeList_[nodeId].rangeBox_.first.second << "] ["
         << nodeList_[nodeId].rangeBox_.second.first << " "
         << nodeList_[nodeId].rangeBox_.second.second << "] "
         << " area="
         << (nodeList_[nodeId].rangeBox_.first.second
             - nodeList_[nodeId].rangeBox_.first.first)
              * (nodeList_[nodeId].rangeBox_.second.second
                 - nodeList_[nodeId].rangeBox_.second.first)
         << " threshold=" << leafMinimumRangeAreaRatio_ * rangeArea_
         << std::endl;
  stream << "[RangeDrivenOctree] Number of cells: "
         << nodeList_[nodeId].cellList_.size() << std::endl;

  return 0;
}

int RangeDrivenOctree::stats(std::ostream &stream) {

  SimplexId leafNumber = 0, nonEmptyLeafNumber = 0, minCellNumber = -1,
            maxCellNumber = -1, storedCellNumber = 0;
  float averageCellNumber = 0;
  SimplexId maxCellId = 0;

  for(size_t i = 0; i < nodeList_.size(); i++) {
    if(!nodeList_[i].childList_.size()) {
      // leaf
      leafNumber++;
      if(nodeList_[i].cellList_.size()) {
        nonEmptyLeafNumber++;
        storedCellNumber += nodeList_[i].cellList_.size();

        averageCellNumber += nodeList_[i].cellList_.size();
        if((minCellNumber == -1)
           || (nodeList_[i].cellList_.size()
               < static_cast<size_t>(minCellNumber)))
          minCellNumber = nodeList_[i].cellList_.size();
        if((maxCellNumber == -1)
           || (nodeList_[i].cellList_.size()
               > static_cast<size_t>(maxCellNumber))) {
          maxCellNumber = nodeList_[i].cellList_.size();
          maxCellId = i;
        }
      }
    }
  }
  averageCellNumber /= nonEmptyLeafNumber;

  stream << "[RangeDrivenOctree] Domain volume: " << domainVolume_ << std::endl;
  stream << "[RangeDrivenOctree] Range area: " << rangeArea_ << std::endl;
  stream << "[RangeDrivenOctree] Leaf number: " << leafNumber << std::endl;
  stream << "[RangeDrivenOctree] Non empty leaf number: " << nonEmptyLeafNumber
         << std::endl;
  stream << "[RangeDrivenOctree] Average cell number: " << averageCellNumber
         << std::endl;
  stream << "[RangeDrivenOctree] Min cell number: " << minCellNumber
         << std::endl;
  stream << "[RangeDrivenOctree] Max cell number: " << maxCellNumber
         << std::endl;
  stream << "[RangeDrivenOctree] Stored cell number: " << storedCellNumber
         << " [input=" << cellNumber_ << "]" << std::endl;

  stream << "[RangeDrivenOctree] Max-cell nodeId: " << maxCellId << std::endl;

  if(debugLevel_ > 5) {
    for(size_t i = 0; i < nodeList_.size(); i++) {
      if(nodeList_[i].cellList_.size())
        statNode(i, stream);
    }
  }

  return 0;
}

bool RangeDrivenOctree::segmentIntersection(
  const std::pair<double, double> &p0,
  const std::pair<double, double> &p1,
  const std::pair<double, double> &q0,
  const std::pair<double, double> &q1) const {

  // NOTE: [q0, q1] is axis aligned
  bool horizontal = true;
  double x = 0, y = 0;

  if(q0.first == q1.first)
    horizontal = false;

  double denP = p1.first - p0.first;
  if(!denP)
    denP = DBL_EPSILON;

  double P = (p1.second - p0.second) / denP;
  if(!P)
    P = DBL_EPSILON;
  double bP = p1.second - P * p1.first;

  if(horizontal) {
    y = q0.second;
    x = (y - bP) / P;
  } else {
    x = q0.first;
    y = P * x + bP;
  }

  // does the intersection lands on the range segment?
  if(horizontal) {
    if((x < fmin(q0.first, q1.first)) || (x > fmax(q0.first, q1.first))) {
      return false;
    }
  } else {
    if((y < fmin(q0.second, q1.second)) || (y > fmax(q0.second, q1.second))) {
      return false;
    }
  }

  // does it lands on the polygon segment's projection?
  if((x < fmin(p0.first, p1.first) || (x > fmax(p0.first, p1.first)))
     && ((y < fmin(p0.second, p1.second)
          || (y > fmax(p0.second, p1.second))))) {
    return false;
  }

  return true;
}
