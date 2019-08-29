/// \ingroup base
/// \class ttk:FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Sept 2016.
///
///\brief TTK processing package that deal with segmentation
// for the merge tree and contour tree
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

#include <algorithm>
#include <iterator>

#include "FTMSegmentation.h"

using namespace std;
using namespace ttk;
using namespace ftm;

// -------
// Segment
// -------

Segment::Segment(SimplexId size) : vertices_(size, nullVertex) {
}

segm_const_it Segment::begin(void) const {
  return vertices_.begin();
}

segm_it Segment::begin(void) {
  return vertices_.begin();
}

void Segment::createFromList(const Scalars *s,
                             list<vector<SimplexId>> &regularList,
                             const bool reverse) {
  auto vectComp = [&](const vector<SimplexId> &a, const vector<SimplexId> &b) {
    return s->isLower(a[0], b[0]);
  };

  SimplexId totalSize = 0;
  for(const auto &vectReg : regularList) {
    totalSize += vectReg.size();
  }
  vertices_.reserve(totalSize);
  // TODO parallel
  regularList.sort(vectComp);
  // TODO parallel
  if(reverse) {
    for(auto it1 = regularList.cbegin(); it1 != regularList.cend(); ++it1) {
      for(auto it2 = it1->crbegin(); it2 != it1->crend(); ++it2) {
        vertices_.emplace_back(*it2);
      }
    }
  } else {
    for(auto it1 = regularList.cbegin(); it1 != regularList.cend(); ++it1) {
      for(auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2) {
        vertices_.emplace_back(*it2);
      }
    }
  }

  regularList.clear();
}

segm_const_it Segment::end(void) const {
  return vertices_.end();
}

segm_it Segment::end(void) {
  return vertices_.end();
}

SimplexId Segment::operator[](const size_t &idx) const {
  return vertices_[idx];
}

SimplexId &Segment::operator[](const size_t &idx) {
  return vertices_[idx];
}

SimplexId Segment::size(void) const {
  return vertices_.size();
}

void Segment::sort(const Scalars *s) {
  // Sort by scalar value
  auto comp = [&](SimplexId a, SimplexId b) { return s->isLower(a, b); };

  std::sort(vertices_.begin(), vertices_.end(), comp);
}

// --------
// Segments
// --------

Segments::Segments() {
}

void Segments::clear(void) {
  segments_.clear();
}

const Segment &Segments::operator[](const size_t &idx) const {
  return segments_[idx];
}

Segment &Segments::operator[](const size_t &idx) {
  return segments_[idx];
}

void Segments::resize(const vector<SimplexId> &sizes) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(segments_.size()) {
    cerr << "Call reserve on an already reserved Segments! " << endl;
  }
#endif

  segments_.reserve(sizes.size());
  for(SimplexId size : sizes) {
    segments_.emplace_back(size);
  }
}

idSegment Segments::size(void) const {
  return segments_.size();
}

void Segments::sortAll(const Scalars *s) {
  const idSegment &nbSegments = size();

  for(idSegment i = 0; i < nbSegments; i++) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i)
#endif
    segments_[i].sort(s);
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

// ----------
// Arc Region
// ----------

ArcRegion::ArcRegion() {
#ifndef TTK_ENABLE_KAMIKAZE
  segmented_ = false;
#endif
}

ArcRegion::ArcRegion(const segm_it &begin, const segm_it &end) : ArcRegion() {
  concat(begin, end);
}

void ArcRegion::concat(const segm_it &begin, const segm_it &end) {
  segmentsIn_.emplace_front(Region{begin, end});
}

void ArcRegion::concat(const ArcRegion &r) {
  for(const auto &reg : r.segmentsIn_) {
    concat(reg.segmentBegin, reg.segmentEnd);
  }
}

void ArcRegion::createSegmentation(const Scalars *s) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(segmentation_.size()) {
    cout << "createSegmentation called on an already segmented region" << endl;
  }
#endif

  SimplexId totalSegmSize = 0;
  vector<segm_const_it> heads, ends;
  for(const auto &region : segmentsIn_) {
    totalSegmSize += distance(region.segmentBegin, region.segmentEnd);
    heads.emplace_back(region.segmentBegin);
    ends.emplace_back(region.segmentEnd);
  }

  segmentation_.clear();
  segmentation_.reserve(
    totalSegmSize); // max size, including discarded vertices

  idSegment nbSegments = heads.size();
  int added = 0;

  while(added != -1) {
    added = -1;
    SimplexId minVert;
    for(idSegment i = 0; i < nbSegments; i++) {
      auto &headIt = heads[i];
      const auto &endIt = ends[i];

      if(headIt == endIt) {
        // end of this area
        heads.erase(heads.begin() + i);
        ends.erase(ends.begin() + i);
        --nbSegments;
        --i;
        continue;
      }

      if(added == -1 || s->isLower(*headIt, minVert)) {
        minVert = *headIt;
        added = i;
      }
    }
    if(added != -1) {
      segmentation_.emplace_back(minVert);
      ++heads[added];
    }
  } // end while

#ifndef TTK_ENABLE_KAMIKAZE
  segmented_ = true;
#endif
}

SimplexId ArcRegion::findBelow(SimplexId v,
                               const Scalars *s,
                               const vector<idCorresp> &vert2treeOther) const {
  // split at v and return remaining vertices

  auto comp = [s](SimplexId a, SimplexId b) { return s->isLower(a, b); };
  SimplexId splitVert = nullVertex;
  const bool chkOther = vert2treeOther.size() > 0;

  for(const auto &reg : segmentsIn_) {
    if(s->isEqLower(*reg.segmentBegin, v)
       && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
      // is v is between beg/end
      // append once
      const auto &oldBeg = reg.segmentBegin;
      auto posV = lower_bound(oldBeg, reg.segmentEnd, v, comp);

      // posV == end() would mean v is not in this range (cf if in for)
      if(posV != oldBeg && *posV != v) {
        --posV;
      }

      if(chkOther) {
        while(posV != oldBeg && vert2treeOther[*posV] == nullCorresp) {
          --posV;
        }
        if(posV == oldBeg && vert2treeOther[*posV] == nullCorresp) {
          continue;
        }
      }

      if(posV != reg.segmentEnd) {
        splitVert = *posV;
      }

    } else if(s->isLower(*(reg.segmentEnd - 1), v)) {
      // reg is completely below v, we give it to remainingRegion
      if(splitVert == nullVertex
         || s->isHigher(*(reg.segmentEnd - 1), splitVert)) {
        auto posV = reg.segmentEnd - 1;
        if(chkOther) {
          const auto &oldBeg = reg.segmentBegin;
          while(posV != oldBeg && vert2treeOther[*posV] == nullCorresp) {
            --posV;
          }
          if(((posV == oldBeg && vert2treeOther[*posV] != nullCorresp)
              || posV != oldBeg)
             && (splitVert == nullVertex || s->isHigher(*posV, splitVert))) {
            splitVert = *posV;
          }
        } else {
          splitVert = *posV;
        }
      }
    }
  }

  return splitVert;
}

bool ArcRegion::merge(const ArcRegion &r) {
  const Region &other = r.segmentsIn_.front();
  Region &self = segmentsIn_.front();

  if(other.segmentBegin == self.segmentEnd) {
    self.segmentEnd = other.segmentEnd;
    return true;
  } else if(other.segmentEnd == self.segmentBegin) {
    self.segmentBegin = other.segmentBegin;
    return true;
  }

  return false;
}

tuple<SimplexId, ArcRegion> ArcRegion::splitBack(SimplexId v,
                                                 const Scalars *s) {
  // split at v and return remaining vertices

  auto comp = [s](SimplexId a, SimplexId b) { return s->isLower(a, b); };
  ArcRegion remainingRegion;
  SimplexId splitVert = nullVertex;

  list<decltype(segmentsIn_)::iterator> willErase;

  for(decltype(segmentsIn_)::iterator it = segmentsIn_.begin();
      it != segmentsIn_.end(); ++it) {
    auto &reg = *it;
    if(s->isEqLower(*reg.segmentBegin, v)
       && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
      // is v is between beg/end
      // append once
      const auto &oldBeg = reg.segmentBegin;
      auto posV = lower_bound(oldBeg, reg.segmentEnd, v, comp);

      // posV == end() would mean v is not in this range (cf if in for)
      if(posV != oldBeg && *posV != v) {
        --posV;
      }

      if(posV != oldBeg) {
        remainingRegion.concat(oldBeg, posV);
      }

      if(posV == reg.segmentEnd) {
        willErase.emplace_back(it);
      } else {
        splitVert = *posV;
        reg.segmentBegin = posV;
      }

    } else if(s->isLower(*(reg.segmentEnd - 1), v)) {
      // reg is completely below v, we give it to remainingRegion
      remainingRegion.concat(reg.segmentBegin, reg.segmentEnd);
      willErase.emplace_back(it);
      if(splitVert == nullVertex
         || s->isHigher(*(reg.segmentEnd - 1), splitVert)) {
        splitVert = *(reg.segmentEnd - 1);
      }
    }
  }

  // remove in this arc segments that have been moved in remaining
  for(auto &tmpIt : willErase) {
    segmentsIn_.erase(tmpIt);
  }

  return make_tuple(splitVert, remainingRegion);
}

tuple<SimplexId, ArcRegion> ArcRegion::splitFront(SimplexId v,
                                                  const Scalars *s) {
  // this function does not create empty region

  auto comp = [s](SimplexId a, SimplexId b) { return s->isLower(a, b); };
  ArcRegion remainingRegion;
  SimplexId splitVert = nullVertex;

  list<decltype(segmentsIn_)::iterator> willErase;

  for(decltype(segmentsIn_)::iterator it = segmentsIn_.begin();
      it != segmentsIn_.end(); ++it) {
    auto &reg = *it;
    if(s->isEqLower(*reg.segmentBegin, v)
       && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
      // is v is between beg/end
      // append once
      const auto &oldEnd = reg.segmentEnd;
      auto posV = lower_bound(reg.segmentBegin, oldEnd, v, comp);

      if(posV != oldEnd) {
        splitVert = *posV;
        remainingRegion.concat(posV, oldEnd);
      }

      if(posV == reg.segmentBegin) {
        willErase.emplace_back(it);
      } else {
        reg.segmentEnd = posV;
      }

    } else if(s->isHigher(*reg.segmentBegin, v)) {
      // reg is completely above v, we give it to remainingRegion
      remainingRegion.concat(reg.segmentBegin, reg.segmentEnd);
      willErase.emplace_back(it);
      if(splitVert == nullVertex || s->isLower(*reg.segmentBegin, splitVert)) {
        // we ignore vertices that does not come frome this arc
        splitVert = *reg.segmentBegin;
      }
    }
  }

  // remove in this arc segments that have been moved in remaining
  for(auto &tmpIt : willErase) {
    segmentsIn_.erase(tmpIt);
  }

  return make_tuple(splitVert, remainingRegion);
}
