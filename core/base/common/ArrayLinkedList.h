/// \ingroup base
/// \class ttk::ArrayLinkedList
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date Mars 2022.
/// \brief This class describes a dynamic size data structure for thread safe
/// computation. It is a linked list of arrays that also stores the current
/// number of elements. Its key feature is that the addition of an element will
/// never cause the moving of the data structure in memory, unlike an
/// std::vector, making the access to an element thread safe even if another
/// thread is adding elements.
///
/// \sa IntegralLines.h %for a usage example.

#include <array>
#include <list>

#pragma once

namespace ttk {
  template <typename datatype, int size>
  class ArrayLinkedList {
  public:
    std::list<std::array<datatype, size>> list_;
    int numberOfElements_;
    // In order to prevent false sharing when creating a
    // std::vector of ArrayLinkedList objects (one element
    // of the std::vector for each thread), it is necessary
    // that one ArrayLinkedList object be bigger than one cache
    // line. Here, we assume a cache line to be 64 bytes.
    std::array<unsigned char, 32> padding_{};
    ArrayLinkedList()
      : list_(std::list<std::array<datatype, size>>({})), numberOfElements_(0) {
    }

    datatype *addArrayElement(datatype element) {
      numberOfElements_ = numberOfElements_ % size;
      if(numberOfElements_ == 0) {
        this->list_.push_back(std::array<datatype, size>({}));
      }
      this->list_.back().at(numberOfElements_) = element;
      this->numberOfElements_++;
      return &(this->list_.back().at(numberOfElements_ - 1));
    }
  };
} // namespace ttk
