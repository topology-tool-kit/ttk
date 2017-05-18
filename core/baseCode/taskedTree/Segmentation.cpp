/*
 * file: Segmentation.cpp
 * description: Segmentation processing package.
 * author: Gueunet Charles
 * date: September 2016
 */

#ifndef withKamikaze
#include <iostream>
#endif

#include <algorithm>
#include <iterator>

#include "Segmentation.h"

using namespace ttk;
using namespace std;

// -------
// Segment
// -------

Segment::Segment(idVertex size) : vertices_(size, nullVertex)
{
}

void Segment::sort(const Scalars* s)
{
   // Sort by scalar value
   auto comp = [&](idVertex a, idVertex b) { return s->isLower(a, b); };

   std::sort(vertices_.begin(), vertices_.end(), comp);
}

void Segment::createFromList(const Scalars* s, list<vector<idVertex>>& regularList, const bool reverse)
{
   auto vectComp = [&](const vector<idVertex>& a, const vector<idVertex>& b) {
      return s->isLower(a[0], b[0]);
   };

   idVertex totalSize = 0;
   for (const auto& vectReg : regularList) {
      totalSize += vectReg.size();
   }
   vertices_.reserve(totalSize);
   // TODO parallel
   regularList.sort(vectComp);
   //TODO parallel
   if (reverse) {
      for (auto it1 = regularList.crbegin(); it1 != regularList.crend(); ++it1) {
         for (auto it2 = it1->crbegin(); it2 != it1->crend(); ++it2) {
            vertices_.emplace_back(*it2);
         }
      }
   } else {
      for (const auto& vectReg : regularList) {
         for (const idVertex v : vectReg) {
            vertices_.emplace_back(v);
         }
      }
   }

   regularList.clear();
}

segm_const_it Segment::begin(void) const
{
   return vertices_.begin();
}

segm_const_it Segment::end(void) const
{
   return vertices_.end();
}

segm_it Segment::begin(void)
{
   return vertices_.begin();
}

segm_it Segment::end(void)
{
   return vertices_.end();
}

idVertex Segment::size(void) const
{
   return vertices_.size();
}

idVertex Segment::operator[](const size_t& idx) const
{
   return vertices_[idx];
}

idVertex& Segment::operator[](const size_t& idx)
{
   return vertices_[idx];
}

// --------
// Segments
// --------

Segments::Segments()
{
}

void Segments::sortAll(const Scalars* s)
{
   const idSegment& nbSegments = size();

   for (idSegment i = 0; i < nbSegments; i++) {
#pragma omp task firstprivate(i)
      segments_[i].sort(s);
   }
#pragma omp taskwait
}

void Segments::resize(const vector<idVertex>& sizes)
{
#ifndef withKamikaze
   if (segments_.size()) {
      cerr << "Call reserve on an already reserved Segments! " << endl;
   }
#endif

   segments_.reserve(sizes.size());
   for (idVertex size : sizes) {
      segments_.emplace_back(size);
   }
}

// vector like

idSegment Segments::size(void) const
{
   return segments_.size();
}

void Segments::clear(void)
{
    segments_.clear();
}

Segment& Segments::operator[](const size_t& idx)
{
   return segments_[idx];
}

const Segment& Segments::operator[](const size_t& idx) const
{
   return segments_[idx];
}

// ---------
// ArcRegion
// ---------

ArcRegion::ArcRegion()
{
#ifndef withKamikaze
   segmented_ = false;
#endif
}

ArcRegion::ArcRegion(const segm_it& begin, const segm_it& end) : ArcRegion()
{
   concat(begin, end);
}

tuple<idVertex, ArcRegion> ArcRegion::splitFront(idVertex v, const Scalars* s)
{
   // this function does not create empty region

   auto comp = [s](idVertex a, idVertex b) { return s->isLower(a, b); };
   ArcRegion remainingRegion;
   idVertex  splitVert = nullVertex;

   list<decltype(segmentsIn_)::iterator> willErase;

   for (decltype(segmentsIn_)::iterator it = segmentsIn_.begin(); it != segmentsIn_.end(); ++it) {
      auto& reg = *it;
      if (s->isEqLower(*reg.segmentBegin, v) && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
         // is v is between beg/end
         // append once
         const auto& oldEnd = reg.segmentEnd;
         auto        posV   = lower_bound(reg.segmentBegin, oldEnd, v, comp);

         if (posV != oldEnd) {
            splitVert = *posV;
            remainingRegion.concat(posV, oldEnd);
         }

         if (posV == reg.segmentBegin) {
            willErase.emplace_back(it);
         } else {
            reg.segmentEnd = posV;
         }

      } else if (s->isHigher(*reg.segmentBegin, v)) {
         // reg is completely above v, we give it to remainingRegion
         remainingRegion.concat(reg.segmentBegin, reg.segmentEnd);
         willErase.emplace_back(it);
         if (splitVert == nullVertex || s->isLower(*reg.segmentBegin, splitVert)) {
            // we ignore vertices that does not come frome this arc
            splitVert = *reg.segmentBegin;
         }
      }
   }

   // remove in this arc segments that have been moved in remaining
   for (auto& tmpIt : willErase) {
      segmentsIn_.erase(tmpIt);
   }

   return make_tuple(splitVert, remainingRegion);
}

tuple<idVertex, ArcRegion> ArcRegion::splitBack(idVertex v, const Scalars* s)
{
   // split at v and return remaining vertices

   auto comp = [s](idVertex a, idVertex b) { return s->isLower(a, b); };
   ArcRegion remainingRegion;
   idVertex  splitVert = nullVertex;

   list<decltype(segmentsIn_)::iterator> willErase;

   for (decltype(segmentsIn_)::iterator it = segmentsIn_.begin(); it != segmentsIn_.end(); ++it) {
      auto& reg = *it;
      if (s->isEqLower(*reg.segmentBegin, v) && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
         // is v is between beg/end
         // append once
         const auto& oldBeg = reg.segmentBegin;
         auto        posV   = lower_bound(oldBeg, reg.segmentEnd, v, comp);

         // posV == end() would mean v is not in this range (cf if in for)
         if (posV != oldBeg && *posV != v) {
            --posV;
         }

         if (posV != oldBeg) {
            remainingRegion.concat(oldBeg, posV);
         }

         if (posV == reg.segmentEnd) {
            willErase.emplace_back(it);
         } else {
            splitVert        = *posV;
            reg.segmentBegin = posV;
         }

      } else if (s->isLower(*(reg.segmentEnd - 1), v)) {
         // reg is completely below v, we give it to remainingRegion
         remainingRegion.concat(reg.segmentBegin, reg.segmentEnd);
         willErase.emplace_back(it);
         if (splitVert == nullVertex || s->isHigher(*(reg.segmentEnd - 1), splitVert)) {
            splitVert = *(reg.segmentEnd - 1);
         }
      }
   }

   // remove in this arc segments that have been moved in remaining
   for (auto& tmpIt : willErase) {
      segmentsIn_.erase(tmpIt);
   }

   return make_tuple(splitVert, remainingRegion);
}

idVertex ArcRegion::findBelow(idVertex v, const Scalars* s,
                              const vector<idCorresp>& vert2treeOther) const
{
   // split at v and return remaining vertices

   auto comp            = [s](idVertex a, idVertex b) { return s->isLower(a, b); };
   idVertex   splitVert = nullVertex;
   const bool chkOther  = vert2treeOther.size() > 0;

   for (const auto& reg : segmentsIn_) {
      if (s->isEqLower(*reg.segmentBegin, v) && s->isEqHigher(*(reg.segmentEnd - 1), v)) {
         // is v is between beg/end
         // append once
         const auto& oldBeg = reg.segmentBegin;
         auto        posV   = lower_bound(oldBeg, reg.segmentEnd, v, comp);

         // posV == end() would mean v is not in this range (cf if in for)
         if (posV != oldBeg && *posV != v) {
            --posV;
         }

         if (chkOther) {
            while (posV != oldBeg && vert2treeOther[*posV] == nullCorresp) {
               --posV;
            }
            if (posV == oldBeg && vert2treeOther[*posV] == nullCorresp) {
               continue;
            }
         }

         if (posV != reg.segmentEnd) {
            splitVert = *posV;
         }

      } else if (s->isLower(*(reg.segmentEnd - 1), v)) {
         // reg is completely below v, we give it to remainingRegion
         if (splitVert == nullVertex || s->isHigher(*(reg.segmentEnd - 1), splitVert)) {
            auto posV = reg.segmentEnd - 1;
            if (chkOther) {
               const auto& oldBeg = reg.segmentBegin;
               while (posV != oldBeg && vert2treeOther[*posV] == nullCorresp) {
                  --posV;
               }
               if (((posV == oldBeg && vert2treeOther[*posV] != nullCorresp) || posV != oldBeg) &&
                   (splitVert == nullVertex || s->isHigher(*posV, splitVert))) {
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

void ArcRegion::concat(const segm_it& begin, const segm_it& end)
{
   segmentsIn_.emplace_front(Region{begin, end});
}

void ArcRegion::concat(const ArcRegion& r)
{
   for (const auto& reg : r.segmentsIn_) {
      concat(reg.segmentBegin, reg.segmentEnd);
   }
}

bool ArcRegion::merge(const ArcRegion& r)
{
   const Region& other = r.segmentsIn_.front();
   Region&       self  = segmentsIn_.front();

   if (other.segmentBegin == self.segmentEnd) {
      self.segmentEnd = other.segmentEnd;
      return true;
   } else if (other.segmentEnd == self.segmentBegin) {
      self.segmentBegin = other.segmentBegin;
      return true;
   }

   return false;
}

void ArcRegion::createSegmentation(const Scalars* s)
{
#ifndef withKamikaze
   if (segmentation_.size()) {
      cout << "createSegmentation called on an already segmented region" << endl;
   }
#endif

   idVertex              totalSegmSize = 0;
   vector<segm_const_it> heads, ends;
   for (const auto& region : segmentsIn_) {
      totalSegmSize += distance(region.segmentBegin, region.segmentEnd);
      heads.emplace_back(region.segmentBegin);
      ends.emplace_back(region.segmentEnd);
   }

   segmentation_.clear();
   segmentation_.reserve(totalSegmSize);  // max size, including discarded vertices

   idSegment nbSegments = heads.size();
   int       added      = 0;

   while (added != -1) {
      added = -1;
      idVertex minVert;
      for (idSegment i = 0; i < nbSegments; i++) {
         auto&       headIt = heads[i];
         const auto& endIt  = ends[i];

         if (headIt == endIt) {
            // end of this area
            heads.erase(heads.begin() + i);
            ends.erase(ends.begin() + i);
            --nbSegments;
            --i;
            continue;
         }

         if (added == -1 || s->isLower(*headIt, minVert)) {
            minVert = *headIt;
            added   = i;
         }
      }
      if (added != -1) {
         segmentation_.emplace_back(minVert);
         ++heads[added];
      }
   }  // end while

#ifndef withKamikaze
   segmented_ = true;
#endif
}
