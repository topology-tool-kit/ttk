/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#include <vector>

namespace ttk::rpd {
  using id_t = int;
  using value_t = double;
  constexpr value_t inf = std::numeric_limits<value_t>::infinity();

  using Simplex = std::vector<id_t>;
  using FiltratedSimplex = std::pair<Simplex, value_t>;
  using PersistencePair= std::pair<FiltratedSimplex, FiltratedSimplex>;
  using Diagram = std::vector<PersistencePair>;
  using MultidimensionalDiagram = std::vector<Diagram>;

  using Edge = std::pair<id_t, id_t>;
  using EdgeSet = std::vector<Edge>;
  using EdgeSetSet = std::vector<EdgeSet>;
}