/// \ingroup base
/// \class ttk::QuadrangulationSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkQuadrangulationSubdivision.cpp % for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>
#include <set>
#include <tuple>

namespace ttk {

  class QuadrangulationSubdivision : public Debug {

  public:
    inline void setSubdivisionLevel(const unsigned int value) {
      subdivisionLevel_ = value;
    }
    inline void setRelaxationIterations(const unsigned int value) {
      relaxationIterations_ = value;
    }
    inline void setLockInputExtrema(const bool value) {
      lockInputExtrema = value;
    }
    inline void setLockAllInputVertices(const bool value) {
      lockAllInputVertices = value;
    }
    inline void setReverseProjection(const bool value) {
      reverseProjection_ = value;
    }
    inline void setShowResError(const bool value) {
      showResError_ = value;
    }
    inline void setHausdorffLevel(const float value) {
      hausdorffLevel_ = value;
    }
    inline void setInputQuads(void *const address, unsigned int size) {
      inputQuads_ = static_cast<Quad *>(address);
      inputQuadNumber_ = size;
    }
    inline void setInputVertices(void *const address, unsigned int size) {
      inputVertices_ = static_cast<Point *>(address);
      inputVertexNumber_ = size;
    }
    inline void setInputVertexIdentifiers(void *const address,
                                          unsigned int size) {
      nearestVertexIdentifier_.resize(size);
      auto inputVertexIdentifiers = static_cast<SimplexId *>(address);
      for(size_t i = 0; i < size; i++) {
        nearestVertexIdentifier_[i] = inputVertexIdentifiers[i];
      }
    }
    inline void setupTriangulation(Triangulation *const triangl) {
      triangulation_ = triangl;
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
        triangulation_->preprocessVertexTriangles();
      }
    }
    int execute();

    inline long long *getQuadBuf() {
      return reinterpret_cast<long long *>(outputQuads_.data());
    }
    inline size_t getQuadNumber() const {
      return outputQuads_.size();
    }
    inline float *getPointsBuf() {
      return reinterpret_cast<float *>(outputPoints_.data());
    }
    inline size_t getPointsNumber() const {
      return outputPoints_.size();
    }

  private:
    // vtkPoint instance with interleaved coordinates (AoS)
    struct Point {
      float x;
      float y;
      float z;
      Point operator+(const Point other) const {
        Point res{};
        res.x = x + other.x;
        res.y = y + other.y;
        res.z = z + other.z;
        return res;
      }
      Point operator*(const float scalar) const {
        Point res{};
        res.x = x * scalar;
        res.y = y * scalar;
        res.z = z * scalar;
        return res;
      }
      Point operator-(Point other) const {
        return *this + other * (-1);
      }
      Point operator/(const float scalar) const {
        return (*this * (1.0F / scalar));
      }
      friend std::ostream &operator<<(std::ostream &stream, const Point &pt) {
        stream << pt.x << " " << pt.y << " " << pt.z;
        return stream;
      }
    };
    // VTK_QUAD representation with vtkIdType
    struct Quad {
      long long n; // number of vertices, 4
      long long i; // index of first vertex
      long long j; // second vertex
      long long k; // third vertex
      long long l; // fourth vertex
    };

    /**
     * @brief Subdivise a quadrangular mesh
     *
     * Add five new points per existing quadrangle:
     * - four at the middle of the four quadrangle edges
     * - one at the quadrangle barycenter
     *
     * The points are generated once, no duplicate is generated when
     * considering two quadrangles sharing an edge.
     *
     * New quadrangles using these new points are then generated
     * out-of-place in an output vector of quadrangles.
     *
     * @return 0 in case of success
     */
    int subdivise();

    /**
     * @brief Project a generated quadrangle vertex into the
     * triangular input mesh
     *
     * Project a subset of the current quadrangular mesh onto the
     * triangular input mesh.
     *
     * @param[in] filtered Set of indices that should not be projected
     * @param[in] lastIter Indicate last projection iteration for
     * post-processing
     * @return 0 in case of success
     */
    int project(const std::set<size_t> &filtered, bool lastIter = false);

    /**
     * @brief Relax every generated point of a quadrangular mesh
     *
     * Take every generated point of the current quadrangular mesh,
     * and move its position to the barycenter of its neighbors.
     *
     * @param[in] filtered Set of indices that should not be projected
     * @return 0 in case of success
     */
    int relax(const std::set<size_t> &filtered);

    /**
     * @brief Store for every quad vertex its neighbors
     *
     * Each quad vertex should be linked to four other vertices. This
     * functions stores into the quadNeighbors_ member this relation.
     *
     * @param[in] quads Quadrangular mesh to find neighbors in
     * @param[in] secondNeighbors Also store secondary neighbors (quad third
     * vertex)
     *
     * @return 0 in case of success
     */
    int getQuadNeighbors(const std::vector<Quad> &quads,
                         std::vector<std::set<size_t>> &neighbors,
                         bool secondNeighbors = false) const;

    /**
     * @brief Compute the normal of the quadrangulation at point a
     *
     * @param[in] a input index of quadrangle vertex
     *
     * @return normal to quad surface at point a
     */
    Point getQuadNormal(size_t a) const;

    /**
     * @brief Compute the projection in the nearest triangle
     *
     * @param[in] a input index of quadrangle vertex
     * @param[in] forceReverseProj Try reverse projection
     *
     * @return (coordinates of projection, nearest vertex id, number
     * of triangles checked for result, projection id)
     */
    std::tuple<Point, SimplexId, size_t, SimplexId>
      findProjection(size_t a, bool forceReverseProj) const;

    /**
     * @brief Find the middle of a quad edge using Dijkstra
     *
     * Minimize the sum of the distance to the two edge vertices, and the
     * distance absolute difference
     *
     * @param[in] a First edge vertex index
     * @param[in] b Second edge vertex index
     *
     * @return TTK identifier of potential edge middle
     */
    SimplexId findEdgeMiddle(size_t a, size_t b) const;

    /**
     * @brief Find a quad barycenter using Dijkstra
     *
     * Minimize the sum of the distance to every vertex of the current quad.
     *
     * @param[in] quadVertices Vector of quad vertices point ids in which to
     * find a barycenter
     *
     * @return TTK identifier of potential barycenter
     */
    SimplexId findQuadBary(const std::vector<size_t> &quadVertices) const;

    /**
     * @brief Find input vertices with more than 4 neighbors
     *
     * @param[out] output Output set of input extraordinary point indices
     *
     * @return 0 in case of success
     */
    int findExtraordinaryVertices(std::set<size_t> &output) const;

    /**
     * @brief Compute statistics on generated quadrangles
     *
     * Computes:
     * - quadrangle area
     * - diagonals ratio
     * - ratio between the shortest and the longest edges
     * - ratio between the smallest and the biggest angles
     */
    void quadStatistics();

    /**
     * @brief Clear buffers
     */
    void clearData();

    // number of vertices in the mesh
    SimplexId vertexNumber_{};

    // wanted number of subdivisions of the input quadrangles
    unsigned int subdivisionLevel_{1};
    // number of relaxation iterations
    unsigned int relaxationIterations_{10};
    // lock input extrema
    bool lockInputExtrema{false};
    // lock all input vertices
    bool lockAllInputVertices{false};
    // projection method
    bool reverseProjection_{false};
    // display result despite error
    bool showResError_{false};
    // Hausdorff warning level
    float hausdorffLevel_{200.F};

    // number of input quadrangles
    unsigned int inputQuadNumber_{};
    // input quadrangles
    Quad *inputQuads_{};

    // number of input points (quad vertices)
    unsigned int inputVertexNumber_{};
    // input quadrangle vertices (3D coordinates)
    Point *inputVertices_{};

    // input triangulation
    Triangulation *triangulation_{};

    // array of output quadrangles
    std::vector<Quad> outputQuads_{};
    // array of output quadrangle vertices
    std::vector<Point> outputPoints_{};
    // array mapping quadrangle neighbors
    std::vector<std::set<size_t>> quadNeighbors_{};
    // array of nearest input vertex TTK identifier
    std::vector<SimplexId> nearestVertexIdentifier_{};
    // holds geodesic distance to every other quad vertex sharing a quad
    std::vector<std::vector<float>> vertexDistance_{};

  public:
    // array of output quadrangle vertex valences
    std::vector<SimplexId> outputValences_{};
    // array of output quadrangle vertex type
    // 0 - input (critical) point
    // 1 - edge middle
    // 2 - quadrangle barycenter
    std::vector<SimplexId> outputVertType_{};
    // array of output vertex subdivision level
    std::vector<SimplexId> outputSubdivision_{};
    // number of triangles checked per quad vertex for the last projection
    std::vector<SimplexId> trianglesChecked_{};
    // last projection success per quad vertex
    // 0 - not projected (critical point)
    // 1 - projection alongside quadrangle normal
    // 2 - projection alongside triangle normal
    // 3 - failed projection
    std::vector<SimplexId> projSucceeded_{};

    // quadrangles statistics
    std::vector<float> quadArea_{};
    std::vector<float> quadDiagsRatio_{};
    std::vector<float> quadEdgesRatio_{};
    std::vector<float> quadAnglesRatio_{};
    std::vector<float> hausdorff_{};
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <QuadrangulationSubdivision.cpp>
