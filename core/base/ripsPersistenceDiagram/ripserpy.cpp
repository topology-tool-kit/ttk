/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// In this file is implemented Ripser.py, an evolution made by Christopher
/// Tralie of Ripser, a C++ code for the computation of Vietoris-Rips
/// persistence barcodes by Ulrich Bauer. The present software has been slightly
/// modified to allow the extraction of persistence values and persistent
/// simplices. Below can be found the licence of the original code.

/// \b Related \b publications \n
/// "Ripser.py: A Lean Persistent Homology Library for Python" \n
/// Christopher Tralie, Nathaniel Saul and Rann Bar-On \n
/// The Open Journal, (2018).
/// "Ripser: efficient computation of Vietoris-Rips persistence barcodes" \n
/// Ulrich Bauer \n
/// Journal of Applied and Computational Topology, (2021).

/*

Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

MIT License

Original Copyright 2015-2018 Ulrich Bauer.
Modifications Copyright 2018 Christopher Tralie


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the author of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.


*/

#include "ripser.h"

using namespace ripser;

void check_overflow(index_t i);
coefficient_t get_modulo(const coefficient_t val, const coefficient_t modulus);
coefficient_t normalize(const coefficient_t n, const coefficient_t modulus);
std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m);

#ifdef USE_ROBINHOOD_HASHMAP

#include "robin_hood.h"

template <class Key, class T, class H, class E>
using hash_map = robin_hood::unordered_map<Key, T, H, E>;
template <class Key>
using hash = robin_hood::hash<Key>;

#else

template <class Key, class T, class H, class E>
using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key>
using hash = std::hash<Key>;

#endif

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif

static const std::string clear_line("\r\033[K");

static const size_t num_coefficient_bits = 8;

#if defined(TTK_ENABLE_RIPSER_128BITS_IDS) \
  && (defined(__GNUC__) || defined(__clang__))
static const index_t max_simplex_index
  = (__int128(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i) {
  if
#ifdef USE_COEFFICIENTS
    (i > max_simplex_index)
#else
    (i < 0)
#endif
    throw std::overflow_error("simplex index " + std::to_string((uint64_t)i)
                              + " in filtration is larger than maximum index");
}
#else
// 1L on windows is ALWAYS 32 bits, when on unix systems is pointer size
static const index_t max_simplex_index
  = (uintptr_t(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i) {
  if
#ifdef USE_COEFFICIENTS
    (i > max_simplex_index)
#else
    (i < 0)
#endif
    throw std::overflow_error("simplex index " + std::to_string((uint64_t)i)
                              + " in filtration is larger than maximum index "
                              + std::to_string(max_simplex_index));
}
#endif

class binomial_coeff_table {
  /* Using flatten matrix */
  std::vector<index_t> B;
  size_t offset;

public:
  binomial_coeff_table(index_t n, index_t k) : B((n + 1) * (k + 1)) {
    offset = k + 1;
    for(index_t i = 0; i <= n; ++i) {
      B[i * offset] = 1;
      for(index_t j = 1; j < std::min(i, k + 1); ++j)
        B[i * offset + j]
          = B[(i - 1) * offset + j - 1] + B[(i - 1) * offset + j];
      if(i <= k)
        B[i * offset + i] = 1;
      check_overflow(B[i * offset + std::min(i >> 1, k)]);
    }
  }

  index_t operator()(index_t n, index_t k) const {
    assert(n < index_t(B.size() / offset) && k < index_t(offset) && n >= k - 1);
    return B[n * offset + k];
  }
};

/* Modulo operator is expensive, using a mask when modulus is equal 2
 * is much less expensive and speed-ups where observed
 */
coefficient_t get_modulo(const coefficient_t val, const coefficient_t modulus) {
  return (modulus == 2) ? val & 1 : val % modulus;
}

coefficient_t normalize(const coefficient_t n, const coefficient_t modulus) {
  return n > modulus / 2 ? n - modulus : n;
}

std::vector<coefficient_t>
  multiplicative_inverse_vector(const coefficient_t m) {
  std::vector<coefficient_t> inverse(m);
  inverse[1] = 1;
  // m = a * (m / a) + m % a
  // Multiplying with inverse(a) * inverse(m % a):
  // 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
  for(coefficient_t a = 2; a < m; ++a)
    inverse[a] = m - (inverse[m % a] * (m / a)) % m;
  return inverse;
}

#ifdef USE_COEFFICIENTS

// https://stackoverflow.com/a/3312896/13339777
#ifdef _MSC_VER
#define PACK(...) __pragma(pack(push, 1)) __VA_ARGS__ __pragma(pack(pop))
#else
#define PACK(...) __attribute__((__packed__)) __VA_ARGS__
#endif

PACK(struct entry_t {
  index_t index : 8 * sizeof(index_t) - num_coefficient_bits;
  index_t coefficient : num_coefficient_bits;
  entry_t(index_t _index, coefficient_t _coefficient)
    : index(_index), coefficient(_coefficient) {
  }
  entry_t(index_t _index) : index(_index), coefficient(0) {
  }
  entry_t() : index(0), coefficient(0) {
  }
});

static_assert(sizeof(entry_t) == sizeof(index_t),
              "size of entry_t is not the same as index_t");

entry_t make_entry(index_t i, coefficient_t c);
index_t get_index(const entry_t &e);
index_t get_coefficient(const entry_t &e);
void set_coefficient(entry_t &e, const coefficient_t c);
std::ostream &operator<<(std::ostream &stream, const entry_t &e);

entry_t make_entry(index_t i, coefficient_t c) {
  return entry_t(i, c);
}
index_t get_index(const entry_t &e) {
  return e.index;
}
index_t get_coefficient(const entry_t &e) {
  return e.coefficient;
}
void set_coefficient(entry_t &e, const coefficient_t c) {
  e.coefficient = c;
}

std::ostream &operator<<(std::ostream &stream, const entry_t &e) {
  stream << get_index(e) << ":" << get_coefficient(e);
  return stream;
}

#else

using entry_t = index_t;
index_t get_index(const entry_t &i);
index_t get_coefficient(const entry_t & /*i*/);
entry_t make_entry(index_t _index, coefficient_t /*_value*/);
void set_coefficient(entry_t & /*e*/, const coefficient_t /*c*/);

index_t get_index(const entry_t &i) {
  return i;
}
index_t get_coefficient(const entry_t & /*i*/) {
  return 1;
}
entry_t make_entry(index_t _index, coefficient_t /*_value*/) {
  return entry_t(_index);
}
void set_coefficient(entry_t & /*e*/, const coefficient_t /*c*/) {
}

#endif

const entry_t &get_entry(const entry_t &e);
const entry_t &get_entry(const entry_t &e) {
  return e;
}

using diameter_index_t = std::pair<value_t, index_t>;
value_t get_diameter(const diameter_index_t &i);
index_t get_index(const diameter_index_t &i);
value_t get_diameter(const diameter_index_t &i) {
  return i.first;
}
index_t get_index(const diameter_index_t &i) {
  return i.second;
}

using index_diameter_t = std::pair<index_t, value_t>;
index_t get_index(const index_diameter_t &i);
value_t get_diameter(const index_diameter_t &i);
index_t get_index(const index_diameter_t &i) {
  return i.first;
}
value_t get_diameter(const index_diameter_t &i) {
  return i.second;
}

struct diameter_entry_t : std::pair<value_t, entry_t> {
  using std::pair<value_t, entry_t>::pair;
  diameter_entry_t() = default;
  diameter_entry_t(value_t _diameter,
                   index_t _index,
                   coefficient_t _coefficient)
    : diameter_entry_t(_diameter, make_entry(_index, _coefficient)) {
  }
  diameter_entry_t(const diameter_index_t &_diameter_index,
                   coefficient_t _coefficient)
    : diameter_entry_t(get_diameter(_diameter_index),
                       make_entry(get_index(_diameter_index), _coefficient)) {
  }
  diameter_entry_t(const index_t &_index) : diameter_entry_t(0, _index, 0) {
  }
};

const entry_t &get_entry(const diameter_entry_t &p);
entry_t &get_entry(diameter_entry_t &p);
index_t get_index(const diameter_entry_t &p);
coefficient_t get_coefficient(const diameter_entry_t &p);
const value_t &get_diameter(const diameter_entry_t &p);
void set_coefficient(diameter_entry_t &p, const coefficient_t c);

const entry_t &get_entry(const diameter_entry_t &p) {
  return p.second;
}
entry_t &get_entry(diameter_entry_t &p) {
  return p.second;
}
index_t get_index(const diameter_entry_t &p) {
  return get_index(get_entry(p));
}
coefficient_t get_coefficient(const diameter_entry_t &p) {
  return get_coefficient(get_entry(p));
}
const value_t &get_diameter(const diameter_entry_t &p) {
  return p.first;
}
void set_coefficient(diameter_entry_t &p, const coefficient_t c) {
  set_coefficient(get_entry(p), c);
}

template <typename Entry>
struct greater_diameter_or_smaller_index {
  bool operator()(const Entry &a, const Entry &b) {
    return (get_diameter(a) > get_diameter(b))
           || ((get_diameter(a) == get_diameter(b))
               && (get_index(a) < get_index(b)));
  }
};

enum compressed_matrix_layout : std::uint8_t {
  LOWER_TRIANGULAR,
  UPPER_TRIANGULAR
};

template <compressed_matrix_layout Layout>
class compressed_distance_matrix {
public:
  std::vector<value_t> distances;
  std::vector<value_t *> rows;

  compressed_distance_matrix(std::vector<value_t> &&_distances)
    : distances(std::move(_distances)),
      rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
    assert(distances.size() == size() * (size() - 1) / 2);
    init_rows();
  }

  template <typename DistanceMatrix>
  compressed_distance_matrix(const DistanceMatrix &mat)
    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
    init_rows();

    for(size_t i = 1; i < size(); ++i)
      for(size_t j = 0; j < i; ++j)
        rows[i][j] = mat(i, j);
  }

  value_t operator()(const index_t i, const index_t j) const;

  size_t size() const {
    return rows.size();
  }
  void init_rows();
};

using compressed_lower_distance_matrix
  = compressed_distance_matrix<LOWER_TRIANGULAR>;
using compressed_upper_distance_matrix
  = compressed_distance_matrix<UPPER_TRIANGULAR>;

template <>
void compressed_lower_distance_matrix::init_rows() {
  value_t *pointer = &distances[0];
  for(size_t i = 1; i < size(); ++i) {
    rows[i] = pointer;
    pointer += i;
  }
}

template <>
void compressed_upper_distance_matrix::init_rows() {
  value_t *pointer = &distances[0] - 1;
  for(size_t i = 0; i < size() - 1; ++i) {
    rows[i] = pointer;
    pointer += size() - i - 2;
  }
}

template <>
value_t compressed_lower_distance_matrix::operator()(const index_t i,
                                                     const index_t j) const {
  return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_upper_distance_matrix::operator()(const index_t i,
                                                     const index_t j) const {
  return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

struct sparse_distance_matrix {
  std::vector<std::vector<index_diameter_t>> neighbors;
  index_t num_edges;

  mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
    neighbor_it;
  mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
    neighbor_end;

  sparse_distance_matrix(
    std::vector<std::vector<index_diameter_t>> &&_neighbors, index_t _num_edges)
    : neighbors(std::move(_neighbors)), num_edges(_num_edges) {
  }

  template <typename DistanceMatrix>
  sparse_distance_matrix(const DistanceMatrix &mat, const value_t threshold)
    : neighbors(mat.size()), num_edges(0) {
    for(size_t i = 0; i < size(); ++i)
      for(size_t j = 0; j < size(); ++j)
        if(i != j && mat(i, j) <= threshold) {
          ++num_edges;
          neighbors[i].emplace_back(j, mat(i, j));
        }
  }

  size_t size() const {
    return neighbors.size();
  }
};

struct euclidean_distance_matrix {
  std::vector<std::vector<value_t>> points;

  euclidean_distance_matrix(std::vector<std::vector<value_t>> &&_points)
    : points(std::move(_points)) {
  }

  value_t operator()(const index_t i, const index_t j) const {
    assert(i < index_t(points.size()));
    assert(j < index_t(points.size()));
    return std::sqrt(std::inner_product(
      points[i].begin(), points[i].end(), points[j].begin(), value_t(),
      std::plus<value_t>(),
      [](value_t u, value_t v) { return (u - v) * (u - v); }));
  }

  size_t size() const {
    return points.size();
  }
};

class union_find {
  std::vector<index_t> parent;
  std::vector<uint8_t> rank;

public:
  union_find(const index_t n) : parent(n), rank(n, 0) {
    for(index_t i = 0; i < n; ++i)
      parent[i] = i;
  }

  index_t find(index_t x) {
    index_t y = x, z;
    while((z = parent[y]) != y)
      y = z;
    while((z = parent[x]) != y) {
      parent[x] = y;
      x = z;
    }
    return z;
  }

  void link(index_t x, index_t y) {
    if((x = find(x)) == (y = find(y)))
      return;
    if(rank[x] > rank[y])
      parent[y] = x;
    else {
      parent[x] = y;
      if(rank[x] == rank[y])
        ++rank[y];
    }
  }
};

template <typename T>
T begin(std::pair<T, T> &p) {
  return p.first;
}
template <typename T>
T end(std::pair<T, T> &p) {
  return p.second;
}

template <typename ValueType>
class compressed_sparse_matrix {
  std::vector<size_t> bounds;
  std::vector<ValueType> entries;

  using iterator = typename std::vector<ValueType>::iterator;
  using iterator_pair = std::pair<iterator, iterator>;

public:
  size_t size() const {
    return bounds.size();
  }

  iterator_pair subrange(const index_t index) {
    return {entries.begin() + (index == 0 ? 0 : bounds[index - 1]),
            entries.begin() + bounds[index]};
  }

  void append_column() {
    bounds.push_back(entries.size());
  }

  void push_back(const ValueType e) {
    assert(0 < size());
    entries.push_back(e);
    ++bounds.back();
  }

  void pop_back() {
    assert(0 < size());
    entries.pop_back();
    --bounds.back();
  }
};

template <typename DistanceMatrix>
class Ripser {
  const DistanceMatrix dist;
  index_t n, dim_max;
  const value_t threshold;
  const float ratio;
  const coefficient_t modulus;
  const binomial_coeff_table binomial_coeff;
  const std::vector<coefficient_t> multiplicative_inverse;
  mutable std::vector<diameter_entry_t> cofacet_entries;

  struct entry_hash {
    std::size_t operator()(const entry_t &e) const {
#if defined(USE_ROBINHOOD_HASHMAP)
      return robin_hood::hash<index_t>()(::get_index(e));
#else
      return std::hash<index_t>()(::get_index(e));
#endif
    }
  };

  struct equal_index {
    bool operator()(const entry_t &e, const entry_t &f) const {
      return ::get_index(e) == ::get_index(f);
    }
  };

  using entry_hash_map = hash_map<entry_t, size_t, entry_hash, equal_index>;

public:
  Ripser(DistanceMatrix &&_dist,
         index_t _dim_max,
         value_t _threshold,
         float _ratio,
         coefficient_t _modulus)
    : dist(std::move(_dist)), n(dist.size()), dim_max(_dim_max),
      threshold(_threshold), ratio(_ratio), modulus(_modulus),
      binomial_coeff(n, dim_max + 2),
      multiplicative_inverse(multiplicative_inverse_vector(_modulus)) {
  }

  index_t
    get_max_vertex(const index_t idx, const index_t k, const index_t n_) const {
    auto top = n_;
    auto bottom = k - 1;
    if(binomial_coeff(top, k) > idx) {
      index_t count = top - bottom;
      index_t step;
      index_t mid;
      while(count > 0) {
        step = count >> 1;
        mid = top - step;
        if(binomial_coeff(mid, k) > idx) {
          top = mid - 1;
          count -= step + 1;
        } else
          count = step;
      }
    }
    return top;
  }

  index_t get_edge_index(const index_t i, const index_t j) const {
    return binomial_coeff(i, 2) + j;
  }

  template <typename OutputIterator>
  OutputIterator get_simplex_vertices(index_t idx,
                                      const index_t dim,
                                      index_t n_,
                                      OutputIterator out) const {
    --n_;
    for(index_t k = dim + 1; k > 0; --k) {
      n_ = get_max_vertex(idx, k, n_);
      *out++ = n_;
      idx -= binomial_coeff(n_, k);
    }
    return out;
  }

  class simplex_coboundary_enumerator;

  void
    assemble_columns_to_reduce(std::vector<diameter_index_t> &simplices,
                               std::vector<diameter_index_t> &columns_to_reduce,
                               entry_hash_map &pivot_column_index,
                               index_t dim) {
#ifdef INDICATE_PROGRESS
    std::cerr << clear_line << "assembling columns" << std::flush;
    std::chrono::steady_clock::time_point next
      = std::chrono::steady_clock::now() + time_step;
#endif
    --dim;
    columns_to_reduce.clear();
    std::vector<diameter_index_t> next_simplices;

    for(diameter_index_t &simplex : simplices) {
      simplex_coboundary_enumerator cofacets(
        diameter_entry_t(simplex, 1), dim, *this);

      while(cofacets.has_next(false)) {
#ifdef INDICATE_PROGRESS
        if(std::chrono::steady_clock::now() > next) {
          std::cerr << clear_line << "assembling " << next_simplices.size()
                    << " columns (processing "
                    << std::distance(&simplices[0], &simplex) << "/"
                    << simplices.size() << " simplices)" << std::flush;
          next = std::chrono::steady_clock::now() + time_step;
        }
#endif
        auto cofacet = cofacets.next();
        if(get_diameter(cofacet) <= threshold) {
          if(dim != dim_max)
            next_simplices.emplace_back(
              get_diameter(cofacet), get_index(cofacet));

          if(pivot_column_index.find(get_entry(cofacet))
             == pivot_column_index.end())
            columns_to_reduce.emplace_back(
              get_diameter(cofacet), get_index(cofacet));
        }
      }
    }

    simplices.swap(next_simplices);

#ifdef INDICATE_PROGRESS
    std::cerr << clear_line << "sorting " << columns_to_reduce.size()
              << " columns" << std::flush;
#endif

    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
              greater_diameter_or_smaller_index<diameter_index_t>());

#ifdef INDICATE_PROGRESS
    std::cerr << clear_line << std::flush;
#endif
  }

  void compute_dim_0_pairs(std::vector<diameter_index_t> &edges,
                           std::vector<diameter_index_t> &columns_to_reduce,
                           std::vector<std::vector<pers_pair_t>> &ph) {
    union_find dset(n);

    edges = get_edges();
    std::sort(edges.rbegin(), edges.rend(),
              greater_diameter_or_smaller_index<diameter_index_t>());
    std::vector<index_t> vertices_of_edge(2);
    for(auto e : edges) {
      get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin());
      index_t u = dset.find(vertices_of_edge[0]),
              v = dset.find(vertices_of_edge[1]);

      if(u != v) {
        dset.link(u, v);
        if(get_diameter(e) != 0)
          ph[0].emplace_back(
            simplex_diam_t{
              {u == dset.find(vertices_of_edge[0]) ? vertices_of_edge[1]
                                                   : vertices_of_edge[0]},
              0.},
            simplex_diam_t{
              {vertices_of_edge[0], vertices_of_edge[1]}, get_diameter(e)});
      } else
        columns_to_reduce.push_back(e);
    }
    std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

    for(index_t i = 0; i < n; ++i) {
      if(dset.find(i) == i)
        ph[0].emplace_back(
          simplex_diam_t{{i}, 0.},
          simplex_diam_t{{-1}, std::numeric_limits<value_t>::infinity()});
    }
  }

  template <typename Column>
  diameter_entry_t pop_pivot(Column &column) {
    diameter_entry_t pivot(-1);
#ifdef USE_COEFFICIENTS
    while(!column.empty()) {
      if(get_coefficient(pivot) == 0)
        pivot = column.top();
      else if(get_index(column.top()) != get_index(pivot))
        return pivot;
      else
        set_coefficient(pivot, get_modulo((get_coefficient(pivot)
                                           + get_coefficient(column.top())),
                                          modulus));
      column.pop();
    }
    return (get_coefficient(pivot) == 0) ? -1 : pivot;
#else
    while(!column.empty()) {
      pivot = column.top();
      column.pop();
      if(column.empty() || get_index(column.top()) != get_index(pivot))
        return pivot;
      column.pop();
    }
    return -1;
#endif
  }

  template <typename Column>
  diameter_entry_t get_pivot(Column &column) {
    diameter_entry_t result = pop_pivot(column);
    if(get_index(result) != -1)
      column.push(result);
    return result;
  }

  template <typename Column>
  diameter_entry_t
    init_coboundary_and_get_pivot(const diameter_entry_t simplex,
                                  Column &working_coboundary,
                                  const index_t &dim,
                                  entry_hash_map &pivot_column_index) {
    bool check_for_emergent_pair = true;
    cofacet_entries.clear();
    simplex_coboundary_enumerator cofacets(simplex, dim, *this);
    while(cofacets.has_next()) {
      diameter_entry_t cofacet = cofacets.next();
      if(get_diameter(cofacet) <= threshold) {
        cofacet_entries.push_back(cofacet);
        if(check_for_emergent_pair
           && (get_diameter(simplex) == get_diameter(cofacet))) {
          if(pivot_column_index.find(get_entry(cofacet))
             == pivot_column_index.end())
            return cofacet;
          check_for_emergent_pair = false;
        }
      }
    }
    for(auto cofacet : cofacet_entries)
      working_coboundary.push(cofacet);
    return get_pivot(working_coboundary);
  }

  template <typename Column>
  void add_simplex_coboundary(const diameter_entry_t simplex,
                              const index_t &dim,
                              Column &working_reduction_column,
                              Column &working_coboundary) {
    working_reduction_column.push(simplex);
    simplex_coboundary_enumerator cofacets(simplex, dim, *this);
    while(cofacets.has_next()) {
      diameter_entry_t cofacet = cofacets.next();
      if(get_diameter(cofacet) <= threshold)
        working_coboundary.push(cofacet);
    }
  }

  template <typename Column>
  void
    add_coboundary(compressed_sparse_matrix<diameter_entry_t> &reduction_matrix,
                   const std::vector<diameter_index_t> &columns_to_reduce,
                   const size_t index_column_to_add,
                   const coefficient_t factor,
                   const size_t &dim,
                   Column &working_reduction_column,
                   Column &working_coboundary) {
    diameter_entry_t column_to_add(
      columns_to_reduce[index_column_to_add], factor);
    add_simplex_coboundary(
      column_to_add, dim, working_reduction_column, working_coboundary);

    for(diameter_entry_t simplex :
        reduction_matrix.subrange(index_column_to_add)) {
      set_coefficient(simplex, get_coefficient(simplex) * factor % modulus);
      add_simplex_coboundary(
        simplex, dim, working_reduction_column, working_coboundary);
    }
  }

  using working_t
    = std::priority_queue<diameter_entry_t,
                          std::vector<diameter_entry_t>,
                          greater_diameter_or_smaller_index<diameter_entry_t>>;

  void compute_pairs(std::vector<diameter_index_t> &columns_to_reduce,
                     entry_hash_map &pivot_column_index,
                     index_t dim,
                     std::vector<std::vector<pers_pair_t>> &ph) {
    compressed_sparse_matrix<diameter_entry_t> reduction_matrix;
    size_t index_column_to_add;

#ifdef INDICATE_PROGRESS
    std::chrono::steady_clock::time_point next
      = std::chrono::steady_clock::now() + time_step;
#endif

    for(size_t index_column_to_reduce = 0;
        index_column_to_reduce < columns_to_reduce.size();
        ++index_column_to_reduce) {
      diameter_entry_t column_to_reduce(
        columns_to_reduce[index_column_to_reduce], 1);
      value_t diameter = get_diameter(column_to_reduce);

      reduction_matrix.append_column();

      working_t working_reduction_column;
      working_t working_coboundary;

      working_reduction_column.push(column_to_reduce);

      diameter_entry_t pivot = init_coboundary_and_get_pivot(
        column_to_reduce, working_coboundary, dim, pivot_column_index);

      while(true) {
#ifdef INDICATE_PROGRESS
        if(std::chrono::steady_clock::now() > next) {
          std::cerr << clear_line << "reducing column "
                    << index_column_to_reduce + 1 << "/"
                    << columns_to_reduce.size() << " (diameter " << diameter
                    << ")" << std::flush;
          next = std::chrono::steady_clock::now() + time_step;
        }
#endif
        if(get_index(pivot) != -1) {
          auto pair = pivot_column_index.find(get_entry(pivot));
          if(pair != pivot_column_index.end()) {
            entry_t other_pivot = pair->first;
            index_column_to_add = pair->second;
            coefficient_t factor
              = modulus
                - get_modulo(
                  get_coefficient(pivot)
                    * multiplicative_inverse[get_coefficient(other_pivot)],
                  modulus);

            add_coboundary(reduction_matrix, columns_to_reduce,
                           index_column_to_add, factor, dim,
                           working_reduction_column, working_coboundary);

            pivot = get_pivot(working_coboundary);
          } else {
            value_t death = get_diameter(pivot);
            if(death > diameter * ratio) {
              std::vector<index_t> vertices_birth(dim + 1),
                vertices_death(dim + 2);
              get_simplex_vertices(
                get_index(column_to_reduce), dim, n, vertices_birth.rbegin());
              get_simplex_vertices(
                get_index(pivot), dim + 1, n, vertices_death.rbegin());
              ph[dim].emplace_back(simplex_diam_t{vertices_birth, diameter},
                                   simplex_diam_t{vertices_death, death});
            }

            pivot_column_index.insert(
              {get_entry(pivot), index_column_to_reduce});

            pop_pivot(working_reduction_column);
            while(true) {
              diameter_entry_t e = pop_pivot(working_reduction_column);

              if(get_index(e) == -1)
                break;
              assert(get_coefficient(e) > 0);
              reduction_matrix.push_back(e);
            }
            break;
          }
        } else {
          std::vector<index_t> vertices_birth(dim + 1);
          get_simplex_vertices(
            get_index(column_to_reduce), dim, n, vertices_birth.rbegin());
          ph[dim].emplace_back(
            simplex_diam_t{vertices_birth, diameter},
            simplex_diam_t{{-1}, std::numeric_limits<value_t>::infinity()});
          break;
        }
      }
    }
#ifdef INDICATE_PROGRESS
    std::cerr << clear_line << std::flush;
#endif
  }

  std::vector<diameter_index_t> get_edges();
  void compute_barcodes(std::vector<std::vector<pers_pair_t>> &ph) {
    std::vector<diameter_index_t> simplices, columns_to_reduce;

    /* prevent cases where dim_max < 0 */
    if(dim_max < 0)
      dim_max = 0;

    compute_dim_0_pairs(simplices, columns_to_reduce, ph);

    for(index_t dim = 1; dim <= dim_max; ++dim) {
      entry_hash_map pivot_column_index;
      pivot_column_index.reserve(columns_to_reduce.size());

      compute_pairs(columns_to_reduce, pivot_column_index, dim, ph);

      if(dim < dim_max)
        assemble_columns_to_reduce(
          simplices, columns_to_reduce, pivot_column_index, dim + 1);
    }
  }
};

template <>
class Ripser<compressed_lower_distance_matrix>::simplex_coboundary_enumerator {
private:
  index_t idx_below, idx_above, v, k;
  std::vector<index_t> vertices;
  const diameter_entry_t simplex;
  const coefficient_t modulus;
  const compressed_lower_distance_matrix &dist;
  const binomial_coeff_table &binomial_coeff;

public:
  simplex_coboundary_enumerator(
    const diameter_entry_t _simplex,
    index_t _dim,
    const Ripser<compressed_lower_distance_matrix> &parent)
    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1),
      k(_dim + 1), vertices(_dim + 1), simplex(_simplex),
      modulus(parent.modulus), dist(parent.dist),
      binomial_coeff(parent.binomial_coeff) {
    parent.get_simplex_vertices(
      get_index(_simplex), _dim, parent.n, vertices.begin());
  }

  bool has_next(bool all_cofacets = true) {
    return (v >= k && (all_cofacets || binomial_coeff(v, k) > idx_below));
  }

  diameter_entry_t next() {
    while((binomial_coeff(v, k) <= idx_below)) {
      idx_below -= binomial_coeff(v, k);
      idx_above += binomial_coeff(v, k + 1);
      --v;
      --k;
      assert(k != -1);
    }
    value_t cofacet_diameter = get_diameter(simplex);
    for(index_t w : vertices)
      cofacet_diameter = std::max(cofacet_diameter, dist(v, w));
    index_t cofacet_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
    coefficient_t cofacet_coefficient
      = (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
    return diameter_entry_t(
      cofacet_diameter, cofacet_index, cofacet_coefficient);
  }
};

template <>
class Ripser<sparse_distance_matrix>::simplex_coboundary_enumerator {
  index_t idx_below, idx_above, k;
  std::vector<index_t> vertices;
  const diameter_entry_t simplex;
  const coefficient_t modulus;
  const sparse_distance_matrix &dist;
  const binomial_coeff_table &binomial_coeff;
  std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
    &neighbor_it;
  std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
    &neighbor_end;
  index_diameter_t neighbor;

public:
  simplex_coboundary_enumerator(const diameter_entry_t _simplex,
                                const index_t _dim,
                                const Ripser<sparse_distance_matrix> &parent)
    : idx_below(get_index(_simplex)), idx_above(0), k(_dim + 1),
      vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus),
      dist(parent.dist), binomial_coeff(parent.binomial_coeff),
      neighbor_it(dist.neighbor_it), neighbor_end(dist.neighbor_end) {
    neighbor_it.clear();
    neighbor_end.clear();

    parent.get_simplex_vertices(idx_below, _dim, parent.n, vertices.rbegin());

    for(auto v : vertices) {
      neighbor_it.push_back(dist.neighbors[v].rbegin());
      neighbor_end.push_back(dist.neighbors[v].rend());
    }
  }

  bool has_next(bool all_cofacets = true) {
    for(auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; it0 != end0;
        ++it0) {
      neighbor = *it0;
      for(size_t idx = 1; idx < neighbor_it.size(); ++idx) {
        auto &it = neighbor_it[idx], end = neighbor_end[idx];
        while(get_index(*it) > get_index(neighbor))
          if(++it == end)
            return false;
        if(get_index(*it) != get_index(neighbor))
          goto continue_outer;
        else
          neighbor = std::max(neighbor, *it);
      }
      while(k > 0 && vertices[k - 1] > get_index(neighbor)) {
        if(!all_cofacets)
          return false;
        idx_below -= binomial_coeff(vertices[k - 1], k);
        idx_above += binomial_coeff(vertices[k - 1], k + 1);
        --k;
      }
      return true;
    continue_outer:;
    }
    return false;
  }

  diameter_entry_t next() {
    ++neighbor_it[0];
    value_t cofacet_diameter
      = std::max(get_diameter(simplex), get_diameter(neighbor));
    index_t cofacet_index
      = idx_above + binomial_coeff(get_index(neighbor), k + 1) + idx_below;
    coefficient_t cofacet_coefficient
      = (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
    return diameter_entry_t(
      cofacet_diameter, cofacet_index, cofacet_coefficient);
  }
};

template <>
std::vector<diameter_index_t>
  Ripser<compressed_lower_distance_matrix>::get_edges() {
  std::vector<diameter_index_t> edges;
  std::vector<index_t> vertices(2);
  for(index_t index = binomial_coeff(n, 2); index-- > 0;) {
    get_simplex_vertices(index, 1, dist.size(), vertices.rbegin());
    value_t length = dist(vertices[0], vertices[1]);
    if(length <= threshold)
      edges.emplace_back(length, index);
  }
  return edges;
}

template <>
std::vector<diameter_index_t> Ripser<sparse_distance_matrix>::get_edges() {
  std::vector<diameter_index_t> edges;
  for(index_t i = 0; i < n; ++i)
    for(auto n_ : dist.neighbors[i]) {
      index_t j = get_index(n_);
      if(i > j)
        edges.emplace_back(get_diameter(n_), get_edge_index(i, j));
    }
  return edges;
}

void ripser::ripser(std::vector<std::vector<value_t>> points,
                    value_t threshold,
                    index_t dim_max,
                    bool distanceMatrix,
                    std::vector<std::vector<pers_pair_t>> &ph) {
  double ratio = 1;
  coefficient_t modulus = 2;

  ph = std::vector<std::vector<pers_pair_t>>(
    dim_max + 1, std::vector<pers_pair_t>(0));

  if(!distanceMatrix) {
    euclidean_distance_matrix eucl_dist(std::move(points));
    sparse_distance_matrix dist(eucl_dist, threshold);
    Ripser<sparse_distance_matrix> ripser(
      std::move(dist), dim_max, threshold, ratio, modulus);
    ripser.compute_barcodes(ph);
  } else {
    compressed_lower_distance_matrix lower_dist(std::move(points[0]));
    sparse_distance_matrix dist(lower_dist, threshold);
    Ripser<sparse_distance_matrix> ripser(
      std::move(dist), dim_max, threshold, ratio, modulus);
    ripser.compute_barcodes(ph);
  }
}
