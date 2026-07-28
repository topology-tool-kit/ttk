/// \ingroup baseCode
/// \class ttk::dcvf::DiscreteVectorField
/// \author Tanner Finken <finkent@arizona.edu>
/// \author Joshua A. Levine <josh@cs.arizona.edu>
/// \date May 2024.
///
/// \brief TTK %discreteVectorField processing package.
///
/// %DiscreteVectorField is a TTK processing package that handles discrete
/// vector field (in the sense of Discrete Morse Theory). The implementation
/// is largely based on Discrete Gradient filter.
///
/// \sa ttk::dcg::DiscreteGradient
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE VIS 2024.
///
/// \sa ttk::TopologicalSkeleton
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <Geometry.h>
#include <SurfaceGeometrySmoother.h>
#include <Triangulation.h>
#include <VisitedMask.h>

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <set>
#include <utility>

using ttk::dcg::Cell;

namespace ttk {

  namespace dcvf {

    enum connectionValue { NULL_CONNECTION = -1, GHOST_CONNECTION = -2 };

    /**
     * @brief Extended Cell structure for processOutwardStars
     */
    struct CellOutExt : Cell {
      explicit CellOutExt(const int dim, const SimplexId id) : Cell{dim, id} {
      }
      explicit CellOutExt(const int dim,
                          const SimplexId id,
                          const std::array<SimplexId, 3> &lowVerts,
                          const std::array<float, 3> &lowVertWeights,
                          const std::array<uint8_t, 3> &faces)
        : Cell{dim, id}, lowVerts_{lowVerts},
          lowVertWeights_{lowVertWeights}, faces_{faces} {
      }

      // ID values for Outward vertices in current Outward star
      // (1 for edges, 2 for triangles, 3 for tetras)
      const std::array<SimplexId, 3> lowVerts_{};
      // (weightValue for) Outward vertices in current Outward star, align with
      // lowVerts_ (1 for edges, 2 for triangles, 3 for tetras)
      const std::array<float, 3> lowVertWeights_{};
      // indices of faces (cells of dimensions dim_ - 1) in Outward star
      // structure, only applicable for triangles (2 edge faces)  and tetras (3
      // triangle faces)
      const std::array<uint8_t, 3> faces_{};
      // if cell has been paired with another in current Outward star
      bool paired_{false};
    };

    /**
     * Type alias for vectors of 3 dimension
     */
    using vectorValue = SurfaceGeometrySmoother::Point;

    /**
     * Compute and manage a discrete vector field of a function on a
     * triangulation. TTK assumes that the input dataset is made of only one
     * connected component.
     */
    class DiscreteVectorField : virtual public Debug {

    public:
      DiscreteVectorField() {
        this->setDebugMsgPrefix("DiscreteVectorField");
#ifdef TTK_ENABLE_MPI
        hasMPISupport_ = true;
#endif
      }

      /**
       * Give whether the simplification takes into account
       * the entire orbit when returning shortcut V-paths.
       */
      void setReverseFullOrbit(bool data) {
        reverseFullOrbit = data;
      }

      /**
       * Compute the initial discrete vector field of the
       * input vector function on the triangulation.
       */
      template <typename dataType, typename triangulationType>
      int buildField(const triangulationType &triangulation);

      /**
       * Compare the values for a given two vertex ids. Return true if a > b.
       * (Meaning from a to b is downhill/ outward so add to OS)
       */
      template <typename dataType, typename triangulationType>
      bool compare(const triangulationType &triangulation,
                   SimplexId vertexA,
                   SimplexId vertexB,
                   float &weightValue) const;

      /**
       * Set the input vector function.
       *
       * The first parameter is a pointer to the vector field buffer
       * (often provided by ttkUtils::GetVoidPointer()), the second
       * one is a timestamp representing the last modification time of
       * the vector field (often provided by vtkObject::GetMTime()).
       */
      inline void setInputVectorField(const void *const data,
                                      const size_t mTime) {
        inputVectorField_ = std::make_pair(data, mTime);
      }
      /**
       * Get the vector values from inputVectorField_ (needs to be set first)
       * of a particular **vertex** id
       */
      template <typename dataType>
      vectorValue getVectorValueAt(SimplexId &vertex) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(inputVectorField_.first == nullptr) {
          this->printErr(
            "Attempting to access vector before setting input vector field");
          return {};
        }
#endif
        vectorValue result;
        const dataType *dataArray
          = reinterpret_cast<const dataType *>(inputVectorField_.first);
        result[0] = static_cast<float>(dataArray[3 * vertex + 0]);
        result[1] = static_cast<float>(dataArray[3 * vertex + 1]);
        result[2] = static_cast<float>(dataArray[3 * vertex + 2]);

        return result;
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

          data->preconditionEdgeTriangles();
          data->preconditionTriangles();
          if(dim >= 2) {
            data->preconditionBoundaryEdges();
          }
          if(dim == 2) {
            data->preconditionCellEdges();
          } else if(dim == 3) {
            data->preconditionBoundaryTriangles();
            data->preconditionVertexTriangles();
            data->preconditionTriangleEdges();
            data->preconditionTriangleStars();
            data->preconditionCellTriangles();
          }
        }
      }

      /**
       * Get the dimensionality of the triangulation.
       */
      int getDimensionality() const;

      /**
       * Get the number of dimensions available for the cells
       * in the triangulation (equal to dimensionality+1).
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
      template <typename dataType, typename triangulationType>
      bool isBoundary(const Cell &cell,
                      const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is a critical point
       * regarding the discrete vector field, false otherwise.
       */
      bool isCellCritical(const int cellDim, const SimplexId cellId) const;
      bool isCellCritical(const Cell &cell) const;

      /**
       * Return the identifier of the cell paired to the cell given
       * by the user in the discrete vector field.
       */
      template <typename triangulationType>
      SimplexId getPairedCell(const Cell &cell,
                              const triangulationType &triangulation,
                              bool isReverse = false) const;

      /**
       * Return the VPath coming from the given cell.
       */
      template <typename dataType, typename triangulationType>
      int getAscendingPath(const Cell &cell,
                           std::vector<Cell> &vpath,
                           const triangulationType &triangulation,
                           const bool stopOnCycle) const;

      /**
       * Return the VPath coming from the given cell for recursively going
       * through cycles in 2D.
       */
      template <typename dataType, typename triangulationType>
      int getAscendingPathRecursive(const Cell &cell,
                                    std::vector<Cell> &vpath,
                                    const triangulationType &triangulation,
                                    std::vector<char> &previousDescPaths,
                                    std::vector<char> &previousAscPaths) const;

      /**
       * Return the VPath terminating at the given cell.
       */
      template <typename dataType, typename triangulationType>
      int getDescendingPath(const Cell &cell,
                            std::vector<Cell> &vpath,
                            const triangulationType &triangulation,
                            const bool stopOnCycle) const;

      /**
       * Return the VPath terminating at the given cell for recursively going
       * through cycles in 2D.
       */
      template <typename dataType, typename triangulationType>
      int getDescendingPathRecursive(const Cell &cell,
                                     std::vector<Cell> &vpath,
                                     const triangulationType &triangulation,
                                     std::vector<char> &previousDescPaths,
                                     std::vector<char> &previousAscPaths) const;

      /**
       * Return the VPath terminating at the given 2-saddle
       *  restricted to the 2-separatrice of the 1-saddle.
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
       * Return the VPath coming from the given 1-saddle
       *  restricted to the 2-separatrice of the 2-saddle.
       */
      template <typename triangulationType>
      void getAscendingPathThroughWall(
        const Cell &saddle1,
        const Cell &saddle2,
        const std::vector<bool> &isVisited,
        std::vector<Cell> *const vpath,
        const triangulationType &triangulation) const;

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
       * Get the vertex id of the most outward flow on
       * the given cell. (This point might not exist for vector field
       * on cells greater than 2 dimensions)
       */
      template <typename dataType, typename triangulationType>
      SimplexId
        getCellGreaterVertex(const Cell c,
                             const triangulationType &triangulation) const;

      /**
       * Get the vertex id of with the most inward flow on
       * the given cell. (This point might not exist for vector field
       * on cells greater than 2 dimensions)
       */
      template <typename dataType, typename triangulationType>
      SimplexId
        getCellLowerVertex(const Cell c,
                           const triangulationType &triangulation) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output std::vectors are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputVectorField_
       */
      template <typename dataType, typename triangulationType>
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
       * The output std::vectors are modified accordingly.
       */
      template <typename dataType, typename triangulationType>
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
       * Set the Cell Connection(discrete vector) to GHOST_CONNECTION
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
       * Build the glyphs representing the discrete vector field.
       */
      template <typename triangulationType>
      int setVectorGlyphs(std::vector<std::array<float, 3>> &points,
                          std::vector<char> &points_pairOrigins,
                          std::vector<char> &cells_pairTypes,
                          std::vector<SimplexId> &cellsIds,
                          std::vector<char> &cellsDimensions,
                          const triangulationType &triangulation) const;

    private:
      /**
       * Type alias for Outward stars of a given cell
       */
      using outwardStarType = std::array<std::vector<CellOutExt>, 4>;

      /**
       * @brief Store the subcomplexes around vertex for which offset
       * at vertex is maximum
       *
       * @param[in] a Vertex Id
       *
       * @return Outward star as 4 sets of cells (0-cells, 1-cells, 2-cells and
       * 3-cells)
       */
      template <typename dataType, typename triangulationType>
      inline void OutwardStar(outwardStarType &os,
                              const SimplexId a,
                              const triangulationType &triangulation) const;

      /**
       * @brief Return the number of unpaired faces of a given cell in
       * a Outward star
       *
       * @param[in] c Input cell
       * @param[in] os Input Outward star
       *
       * @return Number of unpaired faces and a face id
       */
      std::pair<size_t, SimplexId>
        numUnpairedFaces(const CellOutExt &c, const outwardStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTriangle(const CellOutExt &c,
                                 const outwardStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTetra(const CellOutExt &c,
                              const outwardStarType &ls) const;

      /**
       * @brief Pair cells into discrete vector field
       *
       * @param[in] alpha Cell of lower dimension
       * @param[in] beta Cell of higher dimension
       */
      template <typename triangulationType>
      inline void pairCells(CellOutExt &alpha,
                            CellOutExt &beta,
                            const triangulationType &triangulation);

      /**
       * Implements the ProcessOutwardStars algorithm from
       * "Localized Evaluation for Constructing
       * Discrete Vector Fields", T. Finken, J. Tierny,
       * J. A. Levine
       */
      template <typename dataType, typename triangulationType>
      int processOutwardStars(const triangulationType &triangulation);

      /**
       * @brief Initialize/Allocate discrete vector field memory
       */
      void initMemory(const AbstractTriangulation &triangulation);

    public:
      /**
       * Compute the difference of function values of a pair of cells.
       */
      template <typename dataType, typename triangulationType>
      float getPersistence(const std::vector<Cell> &vpath,
                           const triangulationType &triangulation) const;

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
       * Reverse the given alternating(changes dimensions) VPath
       */
      template <typename triangulationType>
      int reverseAlternatingPath(const std::vector<Cell> &vpath,
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
      bool reverseFullOrbit{true};

      // spare storage for discrete vectors internal structure
      std::array<std::vector<SimplexId>, 6> localVectors_{};
      // former cache key (vector field pointer + timestamp)
      std::pair<const void *, size_t> inputVectorField_{};
      // pointer to either cache entry corresponding to inputVectorField_ or
      // localVectors_ (because cache is not implemented)
      std::array<std::vector<SimplexId>, 6> *vectors_{};
    };

  } // namespace dcvf
} // namespace ttk

#include <DiscreteVectorField_Template.h>
