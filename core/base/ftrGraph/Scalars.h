#pragma once

#include <iostream>
#include <memory>
#include <vector>

#ifdef __APPLE__
#include <algorithm>
#include <numeric>
#else
#ifdef _WIN32
#include <algorithm>
#include <numeric>
#else
#ifdef __clang__
#include <algorithm>
#include <numeric>
#else
#include <parallel/algorithm>
#endif
#endif
#endif

#include <Debug.h>
#include "DataTypes.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      struct Vert {
         idVertex   id;
         ScalarType value;
         idVertex   offset;

         inline bool operator<(const Vert<ScalarType>& v) const
         {
            return value < v.value || (value == v.value && offset < v.offset);
         }
      };

      template <typename ScalarType>
      class Scalars : virtual public Debug
      {
        private:
         idVertex size_;

         ScalarType* values_;
         idVertex*   offsets_;

         bool externalOffsets_;

         std::vector<Vert<ScalarType>> vertices_;
         std::vector<idVertex>         mirror_;

        public:
         Scalars()
             : size_(nullVertex),
               values_(nullptr),
               offsets_(nullptr),
               externalOffsets_(false),
               vertices_(),
               mirror_()
         {
         }

         // Heavy, prevent using it
         Scalars(const Scalars& o) = delete;

         virtual ~Scalars()
         {
            if (!externalOffsets_) {
               delete[] offsets_;
            }
         }

         idVertex getSize(void) const
         {
            return size_;
         }

         ScalarType getVal(const idVertex i)
         {
            return values_[i];
         }

         idVertex getSortedVert(const idVertex i)
         {
            return vertices_[i].id;
         }

         idVertex getMirror(const idVertex i)
         {
            return mirror_[i];
         }

         void setSize(const idVertex size)
         {
            size_ = size;
         }

         void setScalars(ScalarType* values)
         {
            values_ = values;
         }

         void setOffsets(idVertex* sos)
         {
            externalOffsets_ = true;
            offsets_         = sos;
         }

         void alloc()
         {
            if (!externalOffsets_) {
               offsets_ = new idVertex[size_];
            }
            vertices_.resize(size_);
            mirror_.resize(size_);
         }

         void init()
         {
            // Create offset array if not given by user
            if (!externalOffsets_) {
#pragma omp parallel for schedule(static, size_ / threadNumber_)
               for (idVertex i = 0; i < size_; i++) {
                  offsets_[i] = i;
               }
            }

               // Copy everything in the main array
#pragma omp parallel for schedule(static, size_ / threadNumber_)
            for (idVertex i = 0; i < size_; i++) {
               vertices_[i].id     = i;
               vertices_[i].value  = values_[i];
               vertices_[i].offset = offsets_[i];
            }
         }

         void sort()
         {
         // Sort the vertices array
#ifdef TTK_ENABLE_OPENMP
#ifdef __clang__
            cout << "Caution, outside GCC, sequential sort" << endl;
            std::sort(vertices_.begin(), vertices_.end());
#else
            __gnu_parallel::sort(vertices_.begin(), vertices_.end());
#endif
#else
            std::sort(vertices_.begin(), vertices_.end());
#endif

            // Fill the mirror array, used for later comparisons
            DebugTimer tt;
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic, 10000)
            for (idVertex i = 0; i < size_; i++) {
               mirror_[vertices_[i].id] = i;
            }
         }

         void removeNaN(void)
         {
            // This section is aimed to prevent un-deterministic results if the data-set
            // have NaN values in it.
            // In this loop, we replace every NaN by a 0 value.
            // Recall: Equals values are distinguished using Simulation of Simplicity in the FTM
            // tree computation Note: Can we detect NaN using vtk ?
            if (std::numeric_limits<ScalarType>::has_quiet_NaN) {
#pragma omp parallel for
               for (idVertex i = 0; i < size_; i++) {
                  if (isnan(values_[i])) {
                     values_[i] = 0;
                  }
               }
            }
         }

         // Need vertices to be sorted : use mirrorVertices.

         bool isLower(const idVertex a, const idVertex b) const
         {
            return mirror_[a] < mirror_[b];
         }
         bool isEqLower(const idVertex a, const idVertex b) const
         {
            return mirror_[a] <= mirror_[b];
         }

         bool isHigher(const idVertex a, const idVertex b) const
         {
            return mirror_[a] > mirror_[b];
         }
         bool isEqHigher(const idVertex a, const idVertex b) const
         {
            return mirror_[a] >= mirror_[b];
         }
      };
   }
}
