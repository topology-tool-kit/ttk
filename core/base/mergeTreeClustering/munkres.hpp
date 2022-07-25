
#pragma once

//    MIT License

//    Copyright (c) 2019 Aaron Michaux

//    Permission is hereby granted, free of charge, to any person obtaining a copy
//    of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//    copies of the Software, and to permit persons to whom the Software is
//    furnished to do so, subject to the following conditions:

//    The above copyright notice and this permission notice shall be included in all
//    copies or substantial portions of the Software.

//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//    SOFTWARE.

// @see Tutorial on Implementation of Munkres' Assignment Algorithm
//      Robert Pilgram, Murray State University

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>

namespace detail
{
template<typename T> struct MunkresData
{
   static constexpr T zero = T(0);
   enum MunkresState : uint8_t { NONE, STAR, PRIME };
   using edge = std::pair<unsigned, unsigned>;

   unsigned n_rows_ = 0;
   unsigned n_cols_ = 0;
   unsigned side_   = 0;
   std::vector<T> data;
   std::vector<MunkresState> marks;
   std::vector<bool> row_mask;
   std::vector<bool> col_mask;

   // ------------------------------------------------------------- Construction
   //
   MunkresData(const unsigned n_rows,
               const unsigned n_cols,
               std::function<T(unsigned r, unsigned c)> edge_cost) noexcept
       : n_rows_(n_rows)
       , n_cols_(n_cols)
       , side_(std::max(n_rows, n_cols))
       , data(side_ * side_)
       , marks(side_ * side_)
       , row_mask(side_)
       , col_mask(side_)
   {
      assert(n_rows_ > 0);
      assert(n_cols_ > 0);

      // Populate weight matrix... keep track of maximum for next step
      T max_val = std::numeric_limits<T>::lowest();
      for(auto r = 0u; r < n_rows; ++r)
         for(auto c = 0u; c < n_cols; ++c) {
            auto val = edge_cost(r, c);
            C(r, c)  = val;
            if(max_val < val) max_val = val;
         }

      // The weight matrix is always square... fill in the empty
      // spots with 'max-val'
      for(auto r = n_rows; r < side(); ++r)
         for(auto c = n_cols; c < side(); ++c) C(r, c) = max_val;
      for(auto c = n_cols; c < side(); ++c)
         for(auto r = n_rows; r < side(); ++r) C(r, c) = max_val;

      // Subtract the minimum from every row and column, which
      // ensures that every row and column has a '0'
      subtract_min_from_all_rows_cols();

      // Set up marks
      std::fill(begin(marks), end(marks), MunkresState::NONE);
      std::fill(begin(row_mask), end(row_mask), false);
      std::fill(begin(col_mask), end(col_mask), false);

      { // ensure evey element is finite
         auto ii = std::find_if(
             begin(data), end(data), [](float x) { return !std::isfinite(x); });
         if(ii != end(data)) {
            std::cerr << "precondition failed: non-finite edge cost" << std::endl;
            assert(false);
         }
      }
   }

   // ---------------------------------------------------------- Getters/Setters
   //
   // Costs
   T& C(int r, int c) noexcept { return data[r * side_ + c]; }
   const T& C(int r, int c) const noexcept { return data[r * side_ + c]; }

   // Marks
   MunkresState& M(int r, int c) noexcept { return marks[r * side_ + c]; }
   const MunkresState& M(int r, int c) const noexcept { return marks[r * side_ + c]; }

   void cover_row(int r) noexcept { row_mask[r] = true; }
   void cover_col(int c) noexcept { col_mask[c] = true; }
   void uncover_row(int r) noexcept { row_mask[r] = false; }
   void uncover_col(int c) noexcept { col_mask[c] = false; }
   bool is_row_covered(int r) const noexcept { return row_mask[r]; }
   bool is_col_covered(int c) const noexcept { return col_mask[c]; }

   unsigned original_cols() const noexcept { return n_cols_; }
   unsigned original_rows() const noexcept { return n_rows_; }
   unsigned side() const noexcept { return side_; }

   // ------------------------------------------ subtract min from all rows cols
   // This prepares the data for the algorithm
   void subtract_min_from_all_rows_cols()
   {
      auto min_val_in_row = [&](unsigned r) -> T {
         auto min_val = C(r, 0);
         for(auto c = 1u; c < side_; ++c)
            if(C(r, c) < min_val) min_val = C(r, c);
         return min_val;
      };

      auto min_val_in_col = [&](unsigned c) -> T {
         auto min_val = C(0, c);
         for(auto r = 1u; r < side_; ++r)
            if(C(r, c) < min_val) min_val = C(r, c);
         return min_val;
      };

      // Minimize each row
      for(auto r = 0u; r < side_; ++r) {
         const auto min_val = min_val_in_row(r);
         for(auto c = 0u; c < side_; ++c) C(r, c) -= min_val;
      }

      // Minimize each col
      for(auto c = 0u; c < side_; ++c) {
         const auto min_val = min_val_in_col(c);
         for(auto r = 0u; r < side_; ++r) C(r, c) -= min_val;
      }
   }

   // ------------------------------------------------------------------- Step 1
   // Iterate over each element...
   // If it's 0, and there's no other zero in row/col, then STAR
   int step1() noexcept
   {
      std::vector<bool> r_mask(side(), false);
      std::vector<bool> c_mask(side(), false);
      for(auto r = 0u; r < side(); ++r) {
         if(r_mask[r]) continue;
         for(auto c = 0u; c < side(); ++c) {
            if(r_mask[r] || c_mask[c]) continue;
            if(C(r, c) == zero) {
               M(r, c)   = STAR;
               r_mask[r] = true;
               c_mask[c] = true;
            }
         }
      }
      return 2;
   }

   // ------------------------------------------------------------------- Step 2
   // Cover each column containing a STAR
   int step2() noexcept
   {
      auto counter = 0u;

      for(auto c = 0u; c < side(); ++c) assert(!is_col_covered(c));

      for(auto r = 0u; r < side(); ++r) {
         for(auto c = 0u; c < side(); ++c) {
            if(is_col_covered(c)) continue;
            if(M(r, c) == STAR) {
               cover_col(c);
               counter++;
            }
         }
      }

      // A complete matching
      if(counter >= side()) return 0;

      return 3;
   }

   // ------------------------------------------------------------------- Step 3
   // Find a uncovered zero and PRIME it.
   // Eventually get to a state where the PRIMEd row contains no STAR zeros
   std::tuple<int, unsigned, unsigned> step3() noexcept
   {
      auto find_uncovered_row_col = [&](unsigned& r, unsigned& c) -> bool {
         for(r = 0; r < side_; ++r)
            if(!is_row_covered(r))
               for(c = 0; c < side_; ++c)
                  if(!is_col_covered(c))
                     if(C(r, c) == zero) return true;
         return false;
      };

      // Find an uncovered zero, and mark it PRIME
      unsigned saved_row = 0, saved_col = 0;
      if(find_uncovered_row_col(saved_row, saved_col))
         M(saved_row, saved_col) = PRIME;
      else
         return std::tuple<int, unsigned, unsigned>{5, saved_row, saved_col}; // all zeros covered

      // If there's a STAR in the PRIMEd row, then:
      for(auto c = 0u; c < side(); ++c) {
         if(M(saved_row, c) == STAR) {
            cover_row(saved_row);             // cover that row
            uncover_col(c);                   // uncover the column
            return std::tuple<int, unsigned, unsigned>{3, saved_row, saved_col}; // and repeat this step
         }
      }

      // There's no STAR in the PRIMEd row, onto "augmenting path"
      return std::tuple<int, unsigned, unsigned>{4, saved_row, saved_col};
   }

   // ------------------------------------------------------------------- Step 4
   // Augmenting path algorithm
   int step4(const unsigned saved_row, const unsigned saved_col) noexcept
   {
      auto find_star_in_col = [&](const unsigned c) -> int {
         for(auto r = 0u; r < side(); ++r)
            if(M(r, c) == STAR) return r;
         return -1; // row not found
      };

      auto find_prime_in_row = [&](const unsigned r) -> int {
         for(auto c = 0u; c < side(); ++c)
            if(M(r, c) == PRIME) return c;
         assert(false); // we should ALWAYS find this column
         return -1;     // col not found
      };

      auto make_path = [&](const edge e0) {
         std::vector<edge> seq;
         seq.reserve(side());
         seq.push_back(e0);
         int r = -1, c = -1;
         while(true) {
            c = seq.back().second;
            r = find_star_in_col(c); // STARed zero in column of PRIMEd back()
            if(r >= 0)
               seq.push_back({r, c}); // Push a STAR edge
            else                      // If it doesn't exist, then the path is done
               break;
            c = find_prime_in_row(r);
            seq.push_back({r, c}); // Push a PRIME edge
         }
         return seq;
      };

      auto augment_path = [&](const std::vector<edge>& seq) {
         // For all edges in sequence:
         //    1. Erase all STARs
         //    2. And convert all PRIMEs to STARs
         for(const auto& e : seq) {
            if(M(e.first, e.second) == STAR)
               M(e.first, e.second) = NONE;
            else if(M(e.first, e.second) == PRIME)
               M(e.first, e.second) = STAR;
         }
      };

      auto erase_primes = [&]() {
         for(auto r = 0u; r < side(); ++r)
            for(auto c = 0u; c < side(); ++c)
               if(M(r, c) == PRIME) M(r, c) = NONE;
      };

      auto clear_covers = [&]() {
         std::fill(begin(row_mask), end(row_mask), false);
         std::fill(begin(col_mask), end(col_mask), false);
      };

      const edge e0{saved_row, saved_col}; // Uncovered primed zero from step3
      auto seq = make_path(e0);
      augment_path(seq);
      erase_primes();
      clear_covers();

      return 2;
   }

   // ------------------------------------------------------------------- Step 5
   // Find the smallest uncovered value, and:
   //   1. Add it to every covered row
   //   2. Subtract it from every uncovered col
   int step5() noexcept
   {
      auto find_min_uncovered_value = [&]() {
         auto minval = std::numeric_limits<T>::max();
         for(auto r = 0u; r < side(); ++r) {
            if(is_row_covered(r)) continue;
            for(auto c = 0u; c < side(); c++) {
               if(is_col_covered(c)) continue;
               if(C(r, c) < minval) minval = C(r, c);
            }
         }
         return minval;
      };
      const auto minval = find_min_uncovered_value();

      for(auto r = 0u; r < side(); ++r) {
         for(auto c = 0u; c < side(); c++) {
            if(is_row_covered(r)) C(r, c) += minval;  // (1) add minval
            if(!is_col_covered(c)) C(r, c) -= minval; // (2) subtract minval
         }
      }

      return 3;
   }

   // -------------------------------------------------------------------- Solve
   //
   std::vector<edge> solve() noexcept
   {
      // The Munkres Algorithm is described as a state machine
      int step           = 1;
      unsigned saved_row = 0, saved_col = 0;
      while(step) {
         switch(step) {
         case 1:
            step = step1(); // => [2]
            break;
         case 2:
            step = step2(); // => [0, 3]
            break;
         case 3:
            std::tie(step, saved_row, saved_col) = step3(); // => [3, 4, 5]
            break;
         case 4:
            step = step4(saved_row, saved_col); // => [2]
            break;
         case 5:
            step = step5(); // => [3]
            break;
         }
      }

      // Collate the results
      std::vector<edge> out;
      out.reserve(side());
      for(auto r = 0u; r < original_rows(); ++r)
         for(auto c = 0u; c < original_cols(); ++c)
            if(M(r, c) == STAR) out.push_back({r, c});
      return out;
   }
};
} // namespace detail

//
// @param n_lhs_verts Number of left-hand-side vertices (in bipartite graph)
// @param n_rhs_verts Number of right-hand-side verices (ibid)
// @param cost Cost between vertices 'l' and 'r'. Use of function to abstract
//             away storage details of input graph.
// @see example.cpp
//
template<typename T>
std::vector<std::pair<unsigned, unsigned>> inline munkres_algorithm(
    const unsigned n_lhs_verts,
    const unsigned n_rhs_verts,
    std::function<T(unsigned l, unsigned r)> cost) noexcept
{
   detail::MunkresData<T> m_dat{n_lhs_verts, n_rhs_verts, cost};
   return m_dat.solve();
}
