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
      idVertex size_;

      ScalarType *values_;
      std::vector<SimplexId> *vOffsets_;
      SimplexId *offsets_;

      bool externalOffsets_;

      std::vector<Vert<ScalarType>> vertices_;
      std::vector<idVertex> mirror_;

    public:
      Scalars()
        : size_(nullVertex), values_(nullptr), vOffsets_(nullptr),
          offsets_(nullptr), externalOffsets_(false), vertices_(), mirror_() {
      }

      // Heavy, prevent using it
      Scalars(const Scalars &o) = delete;

      virtual ~Scalars() {
        if(!externalOffsets_) {
          delete[] offsets_;
        }
      }

      ScalarType *getScalars() {
        return values_;
      }

      std::vector<SimplexId> *getVOffsets() {
        return vOffsets_;
      }

      SimplexId *getOffsets() {
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

      void setOffsets(std::vector<SimplexId> *sos) {
        externalOffsets_ = sos;
        vOffsets_ = sos;
        if(vOffsets_) {
          offsets_ = vOffsets_->data();
        }
      }

      void alloc() {
        if(!externalOffsets_) {
          offsets_ = new SimplexId[size_];
        }
        vertices_.resize(size_);
        mirror_.resize(size_);
      }

      void init() {
        // Create offset array if not given by user
        if(!externalOffsets_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static, size_ / threadNumber_)
#endif
          for(SimplexId i = 0; i < size_; i++) {
            offsets_[i] = i;
          }
        }

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
        DebugTimer tt;
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
