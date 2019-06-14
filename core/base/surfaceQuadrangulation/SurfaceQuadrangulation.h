/// \ingroup base
/// \class ttk::SurfaceQuadrangulation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkSurfaceQuadrangulation.cpp % for a usage example.

#pragma once

#include <set>

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class SurfaceQuadrangulation : public Debug {

  public:
    inline void setCriticalPoints(const unsigned int number,
                                  void *const points,
                                  void *const ids,
                                  void *const cellIds,
                                  void *const type) {
      criticalPointsNumber_ = number;
      criticalPoints_ = static_cast<float *>(points);
      criticalPointsIdentifier_ = static_cast<SimplexId *>(ids);
      criticalPointsCellIds_ = static_cast<SimplexId *>(cellIds);
      criticalPointsType_ = static_cast<unsigned char *>(type);
    }

    inline void setSeparatrices(const unsigned int number,
                                void *const cellIds,
                                void *const cellDims,
                                void *const mask,
                                void *const points) {
      separatriceNumber_ = number;
      sepCellIds_ = static_cast<SimplexId *>(cellIds);
      sepCellDims_ = static_cast<unsigned char *>(cellDims);
      sepMask_ = static_cast<unsigned char *>(mask);
      sepPoints_ = static_cast<float *>(points);
    }

    inline void setSegmentation(unsigned int number, void *address) {
      segmentationNumber_ = number;
      segmentation_ = static_cast<SimplexId *>(address);
    }
    inline void setDualQuadrangulation(const bool input) {
      dualQuadrangulation_ = input;
    }
    inline void setupTriangulation(Triangulation *const triangl) {
      triangulation_ = triangl;
      if(triangulation_ != nullptr) {
        triangulation_->preprocessVertexNeighbors();
        triangulation_->preprocessVertexTriangles();
      }
    }

    int execute();

  private:
    /**
     * @brief Find the middle of the separatrix specified by its bounds
     *
     * @param[in] a Index in separatrices array of separatrix source
     * @param[in] b Index in separatrices array of separatrix destination
     *
     * @return Index of separatrice source
     */
    size_t findSeparatrixMiddle(const size_t a, const size_t b);

    /**
     * @brief Perform the quadrangulation
     *
     * The direct quadrangulation links extrema to saddle points to
     * make quadrangles.
     *
     * @param[out] ndegen number of degenerate quadrangles produced
     * @return 0 in case of success
     */
    int quadrangulate(size_t &ndegen);

    /**
     * @brief Perform the dual quadrangulation
     *
     * The dual quadrangulation uses only extrema and no saddle points
     * to output a coarser quadrangulation.
     *
     * @return 0 in case of success
     */
    int dualQuadrangulate();

    /**
     * @brief Subdivise quadrangulation
     *
     * Find duplicate separatrices coming from the same vertices and
     * generate new quads that try to map tubular topologies.
     *
     * @return 0 in case of success
     */
    int subdivise();

    /**
     * @brief Find the Morse-Smale cell around a given vertex
     *
     * Perform a breadth-first search from a vertex limited by
     * Morse-Smale separatrices.
     *
     * @param[in] src Source vertex to begin the iteration
     * @param[in,out] vertexCells Cell identifier linked to every vertex (or -1)
     * @param[out] cellSeps Separatrices bordering the current cell
     * @param[in] vertexSepMask If a vertex in on a separatrix
     *
     * @return 0
     */
    int detectCells(const SimplexId src,
                    std::vector<SimplexId> &vertexCells,
                    std::vector<std::vector<SimplexId>> &cellSeps,
                    const std::vector<SimplexId> &vertexSepMask);

    /** @brief Find separatrix index from vertices
     *
     * @param[in] src Source index in critical points array
     * @param[in] dst Destination index in critical points array
     *
     * @return separatrix index
     */
    size_t sepFromPoints(const long long src, const long long dst) const;

    Triangulation *triangulation_{};

    // number of critical points from the Morse-Smale complex
    SimplexId criticalPointsNumber_{};
    // critical points 3d coordinates
    float *criticalPoints_{};
    // mapping points id -> cells id
    SimplexId *criticalPointsCellIds_{};
    // mapping point id -> TTK identifier
    SimplexId *criticalPointsIdentifier_{};
    // critical point type: 0 minimum, 1 saddle point, 2 maximum
    unsigned char *criticalPointsType_{};

    // number of separatrices data
    SimplexId separatriceNumber_{};
    // separatrices points cellIds (to be linked to critical points cellIds)
    SimplexId *sepCellIds_{};
    // separatrices mask scalar field (0 for critical points, 1 otherwise)
    unsigned char *sepMask_{};
    // separatrices cell dimension: 0 for vertices, 1 for edges, 2 for triangles
    unsigned char *sepCellDims_{};
    // separatrices points
    float *sepPoints_{};
    // number of vertices in segmentation
    unsigned int segmentationNumber_{};
    // TTK identifiers -> quad for every vertex segmentation
    SimplexId *segmentation_{};
    // if dual quadrangulation
    bool dualQuadrangulation_{false};

    // index of separatrices beginnings in separatrices arrays
    std::vector<size_t> sepBegs_{};
    // index of separatrices endings in separatrices arrays
    std::vector<size_t> sepEnds_{};
    // sub-segmentation of Morse-Smale cells
    std::vector<SimplexId> morseSeg_{};
    // for each cell, the corresponding MorseSmaleManifold index
    std::vector<SimplexId> cellId_{};
    // indices of separatrices that border quads
    std::vector<std::vector<size_t>> quadSeps_{};

  public:
    // array of output polygons
    std::vector<long long> outputCells_{};
    // array of output vertices (generated middles of duplicated separatrices)
    std::vector<float> outputPoints_{};
    // array of output vertices identifiers
    std::vector<SimplexId> outputPointsIds_{};
    // 0: critical points, 1: edge middle, 2: quad barycenter
    std::vector<SimplexId> outputPointsTypes_{};
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <SurfaceQuadrangulation.cpp>
