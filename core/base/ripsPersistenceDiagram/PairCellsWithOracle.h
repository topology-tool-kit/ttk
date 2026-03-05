/// \ingroup base
/// \class ttk::rpd::PairCellsWithOracle
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date February 2025.
///
/// \brief TTK base class that partially executes on a Rips complex the
/// PairCells persistence algorithm where needed, when the persistence pairs are
/// already known.
///
/// This module defines the %PairCellsWithOracle class that takes a point cloud
/// and enables to compute the persistent generators and the persistent cascades
/// of its Rips complex, when the persistence pairs have already been computed

#pragma once

#include <Debug.h>
#include <RipsPersistenceDiagramUtils.h>
#include <boost/functional/hash.hpp>
#include <set>

namespace ttk::rpd {

  class PairCellsWithOracle : virtual public Debug {
  public:
    PairCellsWithOracle(const PointCloud &points,
                        MultidimensionalDiagram const &oracle,
                        bool distanceMatrix = false,
                        bool parallelSort = false);
    PairCellsWithOracle(float *data,
                        int n,
                        int dim,
                        MultidimensionalDiagram const &oracle,
                        bool parallelSort = false);

    void run();

    void getGenerators(std::vector<Generator1> &generators) const;

    void getCascades(std::vector<Cascade> &cascades, EdgeSets3 &critical) const;
    void getCascades(EdgeSets4 &critical) const;

    static void callOracle(const PointCloud &points,
                           MultidimensionalDiagram &oracle,
                           double threshold = inf,
                           bool distanceMatrix = false);

  private:
    const int n_;
    std::vector<double> compressedDM_{};
    const bool parallelSort_;

    MultidimensionalDiagram const &oracle_;
    double bound_{0.};

    std::unordered_map<Edge, id_t, boost::hash<Edge>> edgeToIndex_;
    std::vector<std::set<id_t, std::greater<>>>
      graph_{}; // std::greater to match the reverse co-lexicographic order used
                // in Ripser
    std::vector<FiltratedEdge> edges_{};
    std::vector<FiltratedTriangle> triangles_{};
    std::vector<id_t> edgesIndices_{};
    std::vector<id_t> edgesOrder_{};

    std::vector<id_t> edgesPartner_{};
    std::vector<id_t> trianglesPartner_{};
    std::vector<std::vector<id_t>> boundaries_{};
    int nEdges_{};

    std::vector<std::vector<id_t>> cascadeEdges_{};

    /**
     * distance matrix access function
     */
    double &DM(unsigned i, unsigned j) {
      return compressedDM_[std::max(i, j) * (std::max(i, j) - 1) / 2
                           + std::min(i, j)];
    }

    void initializeWithBound();
    void pairCellsWithOracle();
    void eliminateBoundaryWithOracle(id_t t_id, id_t e_id);

    template <typename EdgeSets>
    void fillRNG(EdgeSets &critical) const {
      for(auto const &[birth, death] : oracle_[0]) {
        if(death.second < inf)
          critical[0].emplace_back(death.first[0], death.first[1]);
      }
      for(auto const &[birth, death] : oracle_[1])
        critical[1].emplace_back(birth.first[0], birth.first[1]);
    }
  };
} // namespace ttk::rpd
