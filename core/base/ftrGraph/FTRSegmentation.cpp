/// \ingroup base
/// \class ttk::ftr::Segment
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-08-02
///
///\brief TTK processing package that deal with segmentation
/// of an arc in the Reeb Graph

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

#include <algorithm>
#include <iterator>

#include "FTRSegmentation.h"

using namespace std;
using namespace ttk;
using namespace ftr;

// -------
// Segment
// -------

Segment::Segment(idVertex size) : vertices_(size, nullVertex) {
}

Segment::Segment() = default;

segm_const_it Segment::begin() const {
  return vertices_.begin();
}

segm_it Segment::begin() {
  return vertices_.begin();
}

segm_const_it Segment::end() const {
  return vertices_.end();
}

segm_it Segment::end() {
  return vertices_.end();
}

idVertex Segment::size() const {
  return vertices_.size();
}

void Segment::reserve(const idVertex size) {
  vertices_.reserve(size);
}

void Segment::emplace_back(const idVertex v) {
  vertices_.emplace_back(v);
}

idVertex Segment::operator[](const size_t &idx) const {
  return vertices_[idx];
}

idVertex &Segment::operator[](const size_t &idx) {
  return vertices_[idx];
}
