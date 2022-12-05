/// \ingroup base
/// \class ttk::ArrayLinkedList
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date Mars 2022.

#include <array>
#include <list>

#pragma once

namespace ttk {
  template <typename datatype, int size>
  class ArrayLinkedList {
  public:
    std::list<std::array<datatype, size>> list;
    int numberOfElement;
    // In order to prevent false sharing when creating a
    // std::vector of ArrayLinkedList objects (one element
    // of the std::vector for each thread), it is necessary
    // that one ArrayLinkedList object be bigger than one cache
    // line. Here, we assume a cache line to be 64 bytes.
    unsigned char padding[32];
    ArrayLinkedList()
      : list(std::list<std::array<datatype, size>>({})), numberOfElement(0) {
    }

    datatype *addArrayElement(datatype element) {
      numberOfElement = numberOfElement % size;
      if(numberOfElement == 0) {
        this->list.push_back(std::array<datatype, size>({}));
      }
      this->list.back().at(numberOfElement) = element;
      this->numberOfElement++;
      return &(this->list.back().at(numberOfElement - 1));
    }
  };
} // namespace ttk
