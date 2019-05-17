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

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class SurfaceQuadrangulation : public Debug {

  public:
    // default constructor
    SurfaceQuadrangulation() = default;
    // default destructor
    ~SurfaceQuadrangulation() override = default;
    // default copy constructor
    SurfaceQuadrangulation(const SurfaceQuadrangulation &) = default;
    // default move constructor
    SurfaceQuadrangulation(SurfaceQuadrangulation &&) = default;
    // default copy assignment operator
    SurfaceQuadrangulation &operator=(const SurfaceQuadrangulation &) = default;
    // default move assignment operator
    SurfaceQuadrangulation &operator=(SurfaceQuadrangulation &&) = default;

    inline void setCriticalPointsNumber(unsigned int value) {
      criticalPointsNumber_ = value;
    }
    inline void setCriticalPointsCellIds(void *address) {
      criticalPointsCellIds_ = static_cast<SimplexId *>(address);
    }
    inline void setCriticalPointsIdentifiers(void *address) {
      criticalPointsIdentifier_ = static_cast<SimplexId *>(address);
    }
    inline void setCriticalPointsType(void *address) {
      criticalPointsType_ = static_cast<unsigned char *>(address);
    }
    inline void setSeparatriceNumber(unsigned int value) {
      separatriceNumber_ = value;
    }
    inline void setSeparatrices(const unsigned int number,
                                void *const cellIds,
                                void *const mask,
                                void *const points) {
      separatriceNumber_ = number;
      sepCellIds_ = static_cast<SimplexId *>(cellIds);
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
    inline void setOutputQuads(std::vector<long long> *cells,
                               std::vector<float> *points) {
      outputCells_ = cells;
      outputPoints_ = points;
    }

    inline void setupTriangulation(Triangulation *const triangl) {
      triangulation_ = triangl;
      if(triangulation_ != nullptr) {
        triangulation_->preprocessVertexNeighbors();
        triangulation_->preprocessVertexTriangles();
      }
    }

    int execute() const;

  private:
    bool hasCommonManifold(const std::vector<size_t> &verts) const;

    /**
     * @brief Find the middle of the separatrix specified by its vertices
     *
     * @param[in] src Index in critical points array of separatrix source
     * @param[in] dst Index in critical points array of separatrix destination
     *
     * @return Index of middle point in triangulation
     */
    SimplexId findSeparatrixMiddle(const size_t src, const size_t dst) const;

    /**
     * @brief Perform the quadrangulation
     *
     * The direct quadrangulation links extrema to saddle points to
     * make quadrangles.
     *
     * @param[in] sepEdges vector of separatrices edges
     * @param[out] ndegen number of degenerate quadrangles produced
     * @return 0 in case of success
     */
    int quadrangulate(
      const std::vector<std::pair<SimplexId, SimplexId>> &sepEdges,
      size_t &ndegen) const;

    /**
     * @brief Perform the dual quadrangulation
     *
     * The dual quadrangulation uses only extrema and no saddle points
     * to output a coarser quadrangulation.
     *
     * @param[in] sepEdges vector of separatrices edges
     * @return 0 in case of success
     */
    int dualQuadrangulate(
      const std::vector<std::pair<SimplexId, SimplexId>> &sepEdges) const;

    Triangulation *triangulation_{};

    // number of critical points from the Morse-Smale complex
    SimplexId criticalPointsNumber_{};
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
    // separatrices points
    float *sepPoints_{};
    // number of vertices in segmentation
    unsigned int segmentationNumber_{};
    // TTK identifiers -> quad for every vertex segmentation
    SimplexId *segmentation_{};
    // if dual quadrangulation
    bool dualQuadrangulation_{false};

    // array of output polygons
    std::vector<long long> *outputCells_{};
    // array of output edges (critical points and more...)
    std::vector<float> *outputPoints_{};
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <SurfaceQuadrangulation.cpp>
