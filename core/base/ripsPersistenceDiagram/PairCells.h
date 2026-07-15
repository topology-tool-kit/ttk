/// \ingroup base
/// \class ttk::rpd::PairCells
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date September 2024.
///
/// \brief TTK base class that executes the PairCells persistence algorithm
/// on a Rips complex.
///
/// This module defines the %PairCells class that takes a point cloud and
/// enables to compute the persistence diagram, the persistent generators and
/// the persistent cascades of its Rips complex.

#pragma once

#include <Debug.h>
#include <RipsPersistenceDiagramUtils.h>

#ifdef TTK_ENABLE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

#include <set>

namespace ttk::rpd {

  class PairCells : virtual public Debug {
  public:
#ifdef TTK_ENABLE_CGAL
    PairCells(const std::vector<CGAL::Epick::Point_2> &points,
              double upperBound = inf,
              bool parallelSort = false,
              bool parallelMatrixConstruction = false);
#endif
    PairCells(const PointCloud &points,
              bool distanceMatrix = false,
              double upperBound = inf,
              bool parallelSort = false,
              bool parallelMatrixConstruction = false);
    PairCells(float *data,
              int n,
              int dim,
              double upperBound = inf,
              bool parallelSort = false,
              bool parallelMatrixConstruction = false);

    void run();

    void getDiagram(MultidimensionalDiagram &diagrams) const;
    void getDiagramAndGenerators(MultidimensionalDiagram &diagrams,
                                 std::vector<Generator1> &generators) const;

    void getCascades(std::vector<Cascade> &cascades, EdgeSets3 &critical) const;
    void getCascades(EdgeSets4 &critical) const;
    void enrichCascades(std::set<Edge> &cascadeSet,
                        EdgeSets4 &critical,
                        std::vector<int> const &globalIndices) const;

  private:
    const int n_;
    std::vector<double> compressedDM_{};
    const double bound_;
    const bool parallelSort_;
    const bool parallelMatrixConstruction_;

    std::vector<FiltratedEdge> edges_{};
    std::vector<FiltratedTriangle> triangles_{};
    std::vector<id_t> edgesIndices_{};
    std::vector<id_t> edgesOrder_{};
    std::vector<id_t> trianglesIndices_{};

    std::vector<id_t> edgesPartner_{};
    std::vector<id_t> trianglesPartner_{};
    std::vector<std::vector<id_t>> boundaries_{};
    int nEdges_{};
    int nPairedEdges_{0};
    int nTriangles_{};

    std::vector<std::vector<id_t>> cascadeEdges_{};

    /**
     * distance matrix access function
     * \pre i < j
     */
    double &DM(unsigned i, unsigned j) {
      return compressedDM_[j * (j - 1) / 2 + i];
    }

    void initialize();
    void initializeWithBound();

    void executeKruskal();
    void symbolicPerturbation(double eps
                              = std::numeric_limits<float>::epsilon());
    void apparentPairs();

    void pairCells();
    id_t eliminateBoundaries(id_t s);
  };
} // namespace ttk::rpd