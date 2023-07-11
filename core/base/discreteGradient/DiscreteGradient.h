/// \ingroup baseCode
/// \class ttk::dcg::DiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Yizhe Wang <wangyizhe3518@gmail.com>
/// \date November 2016.
///
/// \brief TTK %discreteGradient processing package.
///
/// %DiscreteGradient is a TTK processing package that handles discrete gradient
/// (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <VisitedMask.h>

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <set>
#include <utility>

namespace ttk {

  namespace dcg {
    /**
     * Basic concept of cell, so it must be able to identify any cell of any
     * dimension.
     */
    struct Cell {
      explicit Cell() = default;

      explicit Cell(const int dim, const SimplexId id) : dim_{dim}, id_{id} {
      }

      inline bool operator==(const Cell &other) const {
        return std::tie(this->dim_, this->id_)
               == std::tie(other.dim_, other.id_);
      }

      inline std::string to_string() const {
        return '{' + std::to_string(this->dim_) + ' '
               + std::to_string(this->id_) + '}';
      }

      int dim_{-1};
      SimplexId id_{-1};
    };

    enum gradientValue { NULL_GRADIENT = -1, GHOST_GRADIENT = -2 };

    /**
     * @brief Extended Cell structure for processLowerStars
     */
    struct CellExt : Cell {
      explicit CellExt(const int dim, const SimplexId id) : Cell{dim, id} {
      }
      explicit CellExt(const int dim,
                       const SimplexId id,
                       const std::array<SimplexId, 3> &lowVerts,
                       const std::array<uint8_t, 3> &faces)
        : Cell{dim, id}, lowVerts_{lowVerts}, faces_{faces} {
      }

      // (order field value on) lower vertices in current lower star
      // (1 for edges, 2 for triangles, 3 for tetras)
      const std::array<SimplexId, 3> lowVerts_{};
      // indices of faces (cells of dimensions dim_ - 1) in lower star
      // structure, only applicable for triangles (2 edge faces)  and tetras (3
      // triangle faces)
      const std::array<uint8_t, 3> faces_{};
      // if cell has been paired with another in current lower star
      bool paired_{false};
    };

    /**
     * Compute and manage a discrete gradient of a function on a triangulation.
     * TTK assumes that the input dataset is made of only one connected
     * component.
     */
    class DiscreteGradient : virtual public Debug {

    public:
      DiscreteGradient() {
        this->setDebugMsgPrefix("DiscreteGradient");
#ifdef TTK_ENABLE_MPI
        hasMPISupport_ = true;
#endif
      }

      /**
       * Compute the initial gradient field of the input scalar function on the
triangulation.
       */
      template <typename triangulationType>
      int buildGradient(const triangulationType &triangulation,
                        bool bypassCache = false);

      /**
       * Set the input scalar function.
       *
       * The first parameter is a pointer to the scalar field buffer
       * (often provided by ttkUtils::GetVoidPointer()), the second
       * one is a timestamp representing the last modification time of
       * the scalar field (often provided by vtkObject::GetMTime()).
       */
      inline void setInputScalarField(const void *const data,
                                      const size_t mTime) {
        inputScalarField_ = std::make_pair(data, mTime);
      }

      /**
       * Preprocess all the required connectivity requests on the triangulation.
       */
      inline void preconditionTriangulation(AbstractTriangulation *const data) {
        if(data != nullptr) {
          const auto dim{data->getDimensionality()};

          data->preconditionBoundaryVertices();
          data->preconditionVertexNeighbors();
          data->preconditionVertexEdges();
          data->preconditionVertexStars();
          data->preconditionEdges();
          data->preconditionEdgeStars();
          if(dim >= 2) {
            data->preconditionBoundaryEdges();
          }
          if(dim == 2) {
            data->preconditionCellEdges();
          } else if(dim == 3) {
            data->preconditionBoundaryTriangles();
            data->preconditionVertexTriangles();
            data->preconditionEdgeTriangles();
            data->preconditionTriangles();
            data->preconditionTriangleEdges();
            data->preconditionTriangleStars();
            data->preconditionCellTriangles();
          }
        }
      }

      static inline void
        clearCache(const AbstractTriangulation &triangulation) {
        triangulation.gradientCache_.clear();
      }

      /**
       * @brief Use local storage instead of cache
       */
      inline void setLocalGradient() {
        this->gradient_ = &this->localGradient_;
      }

      /**
       * Set the input offset function.
       *
       * @pre For this function to behave correctly in the absence of
       * the VTK wrapper, ttk::preconditionOrderArray() needs to be
       * called to fill the @p data buffer prior to any
       * computation (the VTK wrapper already includes a mechanism to
       * automatically generate such a preconditioned buffer).
       * @see examples/c++/main.cpp for an example use.
       */
      inline void setInputOffsets(const SimplexId *const data) {
        inputOffsets_ = data;
      }

      /**
       * Get the dimensionality of the triangulation.
       */
      int getDimensionality() const;

      /**
       * Get the number of dimensions available for the cells in the
triangulation (equal to dimensionality+1).
       */
      int getNumberOfDimensions() const;

      /**
       * Get the number of cells of the given dimension.
       */
      template <typename triangulationType>
      SimplexId getNumberOfCells(const int dimension,
                                 const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is at boundary, false otherwise.
       */
      template <typename triangulationType>
      bool isBoundary(const Cell &cell,
                      const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is a critical point regarding the
discrete gradient, false otherwise.
       */
      bool isCellCritical(const int cellDim, const SimplexId cellId) const;
      bool isCellCritical(const Cell &cell) const;

      /**
       * Return the identifier of the cell paired to the cell given by the user
in the gradient.
       */
      template <typename triangulationType>
      SimplexId getPairedCell(const Cell &cell,
                              const triangulationType &triangulation,
                              bool isReverse = false) const;

      /**
       * Return the VPath coming from the given cell.
       */
      template <typename triangulationType>
      int getAscendingPath(const Cell &cell,
                           std::vector<Cell> &vpath,
                           const triangulationType &triangulation,
                           const bool enableCycleDetector = false) const;

      /**
       * Return the VPath terminating at the given cell.
       */
      template <typename triangulationType>
      int getDescendingPath(const Cell &cell,
                            std::vector<Cell> &vpath,
                            const triangulationType &triangulation) const;

      /**
       * Return the VPath terminating at the given 2-saddle restricted to the
2-separatrice of the 1-saddle.
       */
      template <typename triangulationType>
      bool getDescendingPathThroughWall(const Cell &saddle2,
                                        const Cell &saddle1,
                                        const std::vector<bool> &isVisited,
                                        std::vector<Cell> *const vpath,
                                        const triangulationType &triangulation,
                                        const bool stopIfMultiConnected = false,
                                        const bool enableCycleDetector
                                        = false) const;

      /**
       * Return the VPath coming from the given 1-saddle restricted to the
2-separatrice of the 2-saddle.
       */
      template <typename triangulationType>
      bool getAscendingPathThroughWall(const Cell &saddle1,
                                       const Cell &saddle2,
                                       const std::vector<bool> &isVisited,
                                       std::vector<Cell> *const vpath,
                                       const triangulationType &triangulation,
                                       const bool stopIfMultiConnected = false,
                                       const bool enableCycleDetector
                                       = false) const;

      /**
       * Return the 2-separatrice terminating at the given 2-saddle.
       */
      template <typename triangulationType>
      int getDescendingWall(const Cell &cell,
                            VisitedMask &mask,
                            const triangulationType &triangulation,
                            std::vector<Cell> *const wall = nullptr,
                            std::vector<SimplexId> *const saddles
                            = nullptr) const;

      /**
       * Return the 2-separatrice coming from the given 1-saddle.
       */
      template <typename triangulationType>
      int getAscendingWall(const Cell &cell,
                           VisitedMask &mask,
                           const triangulationType &triangulation,
                           std::vector<Cell> *const wall = nullptr,
                           std::vector<SimplexId> *const saddles
                           = nullptr) const;

      /**
       * Get the vertex id of with the maximum scalar field value on
       * the given cell.
       */
      template <typename triangulationType>
      SimplexId
        getCellGreaterVertex(const Cell c,
                             const triangulationType &triangulation) const;

      /**
       * Get the vertex id of with the minimum scalar field value on
       * the given cell.
       */
      template <typename triangulationType>
      SimplexId
        getCellLowerVertex(const Cell c,
                           const triangulationType &triangulation) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output vectors are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename triangulationType>
      int setCriticalPoints(
        const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
        std::vector<std::array<float, 3>> &points,
        std::vector<char> &cellDimensions,
        std::vector<SimplexId> &cellIds,
        std::vector<char> &isOnBoundary,
        std::vector<SimplexId> &PLVertexIdentifiers,
        const triangulationType &triangulation) const;

      /**
       * Detect the critical points and build their geometric embedding.
       * The output vectors are modified accordingly.
       */
      template <typename triangulationType>
      int setCriticalPoints(std::vector<std::array<float, 3>> &points,
                            std::vector<char> &cellDimensions,
                            std::vector<SimplexId> &cellIds,
                            std::vector<char> &isOnBoundary,
                            std::vector<SimplexId> &PLVertexIdentifiers,
                            const triangulationType &triangulation) const;

      /**
       * Get the output critical points as a STL vector of cells.
       */
      template <typename triangulationType>
      int getCriticalPoints(
        std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
        const triangulationType &triangulation) const;

#ifdef TTK_ENABLE_MPI
      /**
       * Set the Cell Gradient to GHOST_GRADIENT
       */
      void setCellToGhost(const int cellDim, const SimplexId cellId);

#endif
      /**
       * Compute manifold size for critical extrema
       */
      int setManifoldSize(
        const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
        const SimplexId *const ascendingManifold,
        const SimplexId *const descendingManifold,
        std::vector<SimplexId> &manifoldSize) const;

      /**
       * Build the glyphs representing the discrete gradient vector field.
       */
      template <typename triangulationType>
      int setGradientGlyphs(std::vector<std::array<float, 3>> &points,
                            std::vector<char> &points_pairOrigins,
                            std::vector<char> &cells_pairTypes,
                            std::vector<SimplexId> &cellsIds,
                            std::vector<char> &cellsDimensions,
                            const triangulationType &triangulation) const;

    private:
      /**
       * Type alias for lower stars of a given cell
       */
      using lowerStarType = std::array<std::vector<CellExt>, 4>;

      /**
       * @brief Store the subcomplexes around vertex for which offset
       * at vertex is maximum
       *
       * @param[in] a Vertex Id
       *
       * @return Lower star as 4 sets of cells (0-cells, 1-cells, 2-cells and
       * 3-cells)
       */
      template <typename triangulationType>
      inline void lowerStar(lowerStarType &ls,
                            const SimplexId a,
                            const SimplexId *const offsets,
                            const triangulationType &triangulation) const;

      /**
       * @brief Return the number of unpaired faces of a given cell in
       * a lower star
       *
       * @param[in] c Input cell
       * @param[in] ls Input lower star
       *
       * @return Number of unpaired faces and a face id
       */
      std::pair<size_t, SimplexId>
        numUnpairedFaces(const CellExt &c, const lowerStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTriangle(const CellExt &c,
                                 const lowerStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTetra(const CellExt &c, const lowerStarType &ls) const;

      /**
       * @brief Return the critical type corresponding to given
       * dimension
       */
      CriticalType criticalTypeFromCellDimension(const int dim) const;

      /**
       * @brief Pair cells into discrete gradient field
       *
       * @param[in] alpha Cell of lower dimension
       * @param[in] beta Cell of higher dimension
       */
      template <typename triangulationType>
      inline void pairCells(CellExt &alpha,
                            CellExt &beta,
                            const triangulationType &triangulation);

      /**
       * Implements the ProcessLowerStars algorithm from "Theory and
       * Algorithms for Constructing Discrete Morse Complexes from
       * Grayscale Digital Images", V. Robins, P. J. Wood,
       * A. P. Sheppard
       */
      template <typename triangulationType>
      int processLowerStars(const SimplexId *const offsets,
                            const triangulationType &triangulation);

      /**
       * @brief Initialize/Allocate discrete gradient memory
       */
      void initMemory(const AbstractTriangulation &triangulation);

    public:
      /**
       * Compute the difference of function values of a pair of cells.
       */
      template <typename dataType, typename triangulationType>
      dataType getPersistence(const Cell &up,
                              const Cell &down,
                              const dataType *const scalars,
                              const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is a minimum regarding the discrete
gradient, false otherwise.
       */
      bool isMinimum(const Cell &cell) const;

      /**
       * Return true if the given cell is a 1-saddle regarding the discrete
gradient, false otherwise.
       */
      bool isSaddle1(const Cell &cell) const;

      /**
       * Return true if the given cell is a 2-saddle regarding the discrete
gradient, false otherwise.
       */
      bool isSaddle2(const Cell &cell) const;

      /**
       * Return true if the given cell is a maximum regarding the discrete
gradient, false otherwise.
       */
      bool isMaximum(const Cell &cell) const;

      /**
       * Reverse the given ascending VPath.
       */
      template <typename triangulationType>
      int reverseAscendingPath(const std::vector<Cell> &vpath,
                               const triangulationType &triangulation) const;

      /**
       * Reverse the given descending VPath.
       */
      template <typename triangulationType>
      int reverseDescendingPath(const std::vector<Cell> &vpath,
                                const triangulationType &triangulation) const;

      /**
       * Reverse the given ascending VPath restricted on a 2-separatrice.
       */
      template <typename triangulationType>
      int reverseAscendingPathOnWall(
        const std::vector<Cell> &vpath,
        const triangulationType &triangulation) const;

      /**
       * Reverse the given descending VPath restricted on a 2-separatrice.
       */
      template <typename triangulationType>
      int reverseDescendingPathOnWall(
        const std::vector<Cell> &vpath,
        const triangulationType &triangulation) const;

    protected:
      int dimensionality_{-1};
      SimplexId numberOfVertices_{};

      // spare storage (bypass cache) for gradient internal structure
      AbstractTriangulation::gradientType localGradient_{};
      // cache key (scalar field pointer + timestamp)
      AbstractTriangulation::gradientKeyType inputScalarField_{};
      // pointer to either cache entry corresponding to inputScalarField_ or
      // localGradient_ (if cache is bypassed)
      AbstractTriangulation::gradientType *gradient_{};
      const SimplexId *inputOffsets_{};
    };

  } // namespace dcg
} // namespace ttk

#include <DiscreteGradient_Template.h>
