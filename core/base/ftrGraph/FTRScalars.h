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
    };

    template <typename ScalarType>
    class Scalars : virtual public Debug {
    private:
      idVertex size_{nullVertex};

      ScalarType *values_{};
      const SimplexId *offsets_{};

      std::vector<Vert<ScalarType>> vertices_{};

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
        return offsets_[i];
      }

      void setSize(const idVertex size) {
        size_ = size;
      }

      void setScalars(ScalarType *values) {
        values_ = values;
      }

      /**
       * @pre For this function to behave correctly in the absence of
       * the VTK wrapper, ttk::preconditionOrderArray() needs to be
       * called to fill the @p sos buffer prior to any
       * computation (the VTK wrapper already includes a mecanism to
       * automatically generate such a preconditioned buffer).
       * @see examples/c++/main.cpp for an example use.
       */
      void setOffsets(const SimplexId *const sos) {
        offsets_ = sos;
      }

      void alloc() {
        vertices_.resize(size_);
      }

      void init() {
        // Copy everything in the main array
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static, size_ / threadNumber_)
#endif
        for(idVertex i = 0; i < size_; i++) {
          // ids and value sorted by offset value
          // vertices_[0] -> global minimum
          // vertices_[size_ - 1] -> global maximum
          vertices_[offsets_[i]].id = i;
          vertices_[offsets_[i]].value = values_[i];
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
        return offsets_[a] < offsets_[b];
      }
      bool isHigher(const idVertex a, const idVertex b) const {
        return offsets_[a] > offsets_[b];
      }
    };
  } // namespace ftr
} // namespace ttk
