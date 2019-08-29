/*
 * file: Segmentation.cpp
 * description: Segmentation processing package.
 * author: Gueunet Charles
 * date: September 2016
 */

#include "DeprecatedSegmentation.h"

#include <algorithm>

using namespace std;
using namespace ttk;
using namespace cf;

// -------
// Segment
// -------

Segment::Segment(const bool order) : ascendingOrder_(order) {
}

bool Segment::isAscending(void) const {
  return ascendingOrder_;
}

void Segment::sort(const Scalars *s) {
  // sort with removal of disabled

  auto ascendingSort
    = [&](const vertex &a, const vertex &b) { return s->isLower(a.id, b.id); };

  auto descendingSort
    = [&](const vertex &a, const vertex &b) { return s->isHigher(a.id, b.id); };

  if(ascendingOrder_) {
    std::sort(vertices_.begin(), vertices_.end(), ascendingSort);
  } else {
    std::sort(vertices_.begin(), vertices_.end(), descendingSort);
  }
}

void Segment::emplace_back(const SimplexId &v) {
  vertices_.emplace_back(vertex{v, nullSuperArc});
}

void Segment::clear(void) {
  vertices_.clear();
  vertices_.shrink_to_fit();
}

SimplexId &Segment::operator[](size_t idx) {
  return vertices_[idx].id;
}

const SimplexId &Segment::operator[](size_t idx) const {
  return vertices_[idx].id;
}

sorted_iterator Segment::sbegin() {
  if(ascendingOrder_) {
    return sorted_iterator(vertices_.begin());
  }

  return sorted_iterator(vertices_.rbegin());
}

sorted_iterator Segment::send() {
  if(ascendingOrder_) {
    return sorted_iterator(vertices_.end());
  }

  return sorted_iterator(vertices_.rend());
}

// --------
// Segments
// --------

idSegment Segments::size(void) const {
  return segments_.size();
}

void Segments::clear(void) {
  segments_.clear();
  segments_.shrink_to_fit();
}

Segment &Segments::operator[](size_t idx) {
  return segments_[idx];
}

const Segment &Segments::operator[](size_t idx) const {
  return segments_[idx];
}

// ---------
// ArcRegion
// ---------

ArcRegion::ArcRegion(const segmentIterator &s) {
  segmentsIn_.emplace_front(Region{s, s, bool()});
}

void ArcRegion::addSegment(const segmentIterator &begin,
                           const segmentIterator &end) {
  segmentsIn_.emplace_front(Region{begin, end, bool()});
}

void ArcRegion::createSegmentation(const idSuperArc &thisArc) {
  SimplexId totalSegmSize = 0;
  vector<sorted_iterator> heads, ends;
  for(const auto &region : segmentsIn_) {
    totalSegmSize += distance(region.segmentBegin, region.segmentEnd);
    heads.emplace_back(region.sbegin());
    heads.emplace_back(region.send());
  }

  segmentation_.clear();
  segmentation_.reserve(
    totalSegmSize); // max size, including discarded vertices

  const auto &nbSegments = heads.size();
  bool added = true;

  while(added) {
    added = false;
    for(unsigned i = 0; i < nbSegments; i++) {
      auto &&head = heads[i];
      auto &&end = ends[i];
      // find next potential vertex
      while(head != end && head->ctArc != thisArc)
        ++head;

      if(!added || 1) {
      } // isLower than actual
      throw "Not finished yet";
    }
  }
}
