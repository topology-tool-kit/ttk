#pragma once

#include "FTRCommon.h"
#include "FTRDataTypes.h"

#include <iostream>
#include <memory>
#include <vector>

#include <Debug.h>

namespace ttk {
  namespace ftr {
    template <typename ScalarType>
    struct Vert {
      idVertex id;
      ScalarType value;
      idVertex offset;

      inline bool operator<(const Vert<ScalarType> &v) const {
        return std::tie(value, offset) < std::tie(v.value, v.offset);
      }
    };

    template <typename ScalarType>
    class Scalars : virtual public Debug {
    private:
      idVertex size_{nullVertex};

      ScalarType *values_{};
      const SimplexId *offsets_{};

      std::vector<Vert<ScalarType>> vertices_{};
      std::vector<idVertex> mirror_{};

    public:
      Scalars() {
      }

      // Heavy, prevent using it
      Scalars(const Scalars &o) = delete;

      ScalarType *getScalars() {
        return values_;
      }

      const SimplexId *getOffsets() const {
        return offsets_;
      }

      idVertex getSize(void) const {
        return size_;
      }

      ScalarType getVal(const idVertex i) const {
        return values_[i];
      }

      idVertex getSortedVert(const idVertex i) const {
        return vertices_[i].id;
      }

      idVertex getMirror(const idVertex i) const {
        return mirror_[i];
      }

      void setSize(const idVertex size) {
        size_ = size;
      }

      void setScalars(ScalarType *values) {
        values_ = values;
      }

      void setOffsets(const SimplexId *const sos) {
        offsets_ = sos;
      }

      void alloc() {
        vertices_.resize(size_);
        mirror_.resize(size_);
      }

      void init() {
        // Copy everything in the main array
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static, size_ / threadNumber_)
#endif
        for(idVertex i = 0; i < size_; i++) {
          vertices_[i].id = i;
          vertices_[i].value = values_[i];
          vertices_[i].offset = offsets_[i];
        }
      }

      void sort() {
        // Sort the vertices array
        ::ttk::ftr::parallel_sort<decltype(vertices_.begin())>(
          vertices_.begin(), vertices_.end());

        // Fill the mirror array, used for later comparisons
        Timer tt;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  schedule(static, size_ / threadNumber_)
#endif
        for(idVertex i = 0; i < size_; i++) {
          mirror_[vertices_[i].id] = i;
        }
      }

      void removeNaN(void) {
        // This section is aimed to prevent un-deterministic results if the
        // data-set have NaN values in it. In this loop, we replace every NaN by
        // a 0 value. Recall: Equals values are distinguished using Simulation
        // of Simplicity in the FTM tree computation Note: Can we detect NaN
        // using vtk ?
        if(std::numeric_limits<ScalarType>::has_quiet_NaN) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  schedule(static, size_ / threadNumber_)
#endif
          for(idVertex i = 0; i < size_; i++) {
            if(::std::isnan((double)values_[i])) {
              values_[i] = 0;
            }
          }
        }
      }

      // Need vertices to be sorted : use mirrorVertices.

      bool isLower(const idVertex a, const idVertex b) const {
        return mirror_[a] < mirror_[b];
      }
      bool isEqLower(const idVertex a, const idVertex b) const {
        return mirror_[a] <= mirror_[b];
      }

      bool isHigher(const idVertex a, const idVertex b) const {
        return mirror_[a] > mirror_[b];
      }
      bool isEqHigher(const idVertex a, const idVertex b) const {
        return mirror_[a] >= mirror_[b];
      }
    };
  } // namespace ftr
} // namespace ttk
