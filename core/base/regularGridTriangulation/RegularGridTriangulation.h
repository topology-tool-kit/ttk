/// \ingroup base
/// \class ttk::RegularGridTriangulation
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date June 2023.
///
/// \brief RegularGridTriangulation is an abstract subclass of
/// ttk::AbstractTriangulation that exposes a common API for triangulation
/// on regular grids. This class is meant to be implemented in
/// ttk::ImplicitTriangulation and ttk::PeriodicImplicitTriangulation.
/// \sa Triangulation

#pragma once

#include <array>

// base code includes
#include <AbstractTriangulation.h>

#ifdef TTK_ENABLE_MPI
#include <memory>
#endif // TTK_ENABLE_MPI

namespace ttk {

  class ImplicitTriangulation;
  class PeriodicImplicitTriangulation;

  class RegularGridTriangulation : public AbstractTriangulation {
    friend class ttk::ImplicitTriangulation;
    friend class ttk::PeriodicImplicitTriangulation;

  public:
    RegularGridTriangulation();
    ~RegularGridTriangulation() override = default;

    virtual int setInputGrid(const float &xOrigin,
                             const float &yOrigin,
                             const float &zOrigin,
                             const float &xSpacing,
                             const float &ySpacing,
                             const float &zSpacing,
                             const SimplexId &xDim,
                             const SimplexId &yDim,
                             const SimplexId &zDim)
      = 0;

    int preconditionDistributedVertices() override;
    // offset coordinates of the local grid inside the metaGrid_
    std::array<SimplexId, 3> localGridOffset_{};
    // hold the neighboring ranks vertex bounding boxes (metaGrid_ coordinates)
    std::vector<std::array<SimplexId, 6>> neighborVertexBBoxes_{};
    // hold the neighboring ranks cells bounding boxes (metaGrid_ coordinates)
    std::vector<std::array<SimplexId, 6>> neighborCellBBoxes_{};
    std::shared_ptr<RegularGridTriangulation> metaGrid_;

    virtual void createMetaGrid(const double *const bounds) = 0;

    RegularGridTriangulation(const RegularGridTriangulation &) = default;
    RegularGridTriangulation(RegularGridTriangulation &&) = default;
    RegularGridTriangulation &operator=(const RegularGridTriangulation &)
      = default;
    RegularGridTriangulation &operator=(RegularGridTriangulation &&) = default;

    SimplexId getVertexGlobalIdInternal(const SimplexId lvid) const override;
    SimplexId getVertexLocalIdInternal(const SimplexId gvid) const override;

    SimplexId getCellGlobalIdInternal(const SimplexId lcid) const override;
    SimplexId getCellLocalIdInternal(const SimplexId gcid) const override;

    SimplexId getEdgeGlobalIdInternal(const SimplexId leid) const override;
    SimplexId getEdgeLocalIdInternal(const SimplexId geid) const override;

    SimplexId getTriangleGlobalIdInternal(const SimplexId ltid) const override;
    SimplexId getTriangleLocalIdInternal(const SimplexId gtid) const override;

    int getVertexRankInternal(const SimplexId lvid) const override;
    int getCellRankInternal(const SimplexId lcid) const override;

  protected:
    std::array<SimplexId, 3> dimensions_; // dimensions
    int dimensionality_;

    virtual void vertexToPosition2d(const SimplexId vertex,
                                    SimplexId p[2]) const = 0;
    virtual void vertexToPosition(const SimplexId vertex,
                                  SimplexId p[3]) const = 0;
    virtual void triangleToPosition2d(const SimplexId triangle,
                                      SimplexId p[2]) const = 0;
    virtual void triangleToPosition(const SimplexId triangle,
                                    const int k,
                                    SimplexId p[3]) const = 0;
    virtual void tetrahedronToPosition(const SimplexId tetrahedron,
                                       SimplexId p[3]) const = 0;

    SimplexId findEdgeFromVertices(const SimplexId v0,
                                   const SimplexId v1) const;
    SimplexId findTriangleFromVertices(std::array<SimplexId, 3> &verts) const;

    virtual std::array<SimplexId, 3>
      getVertGlobalCoords(const SimplexId lvid) const = 0;
    virtual std::array<SimplexId, 3>
      getVertLocalCoords(const SimplexId gvid) const = 0;
  };
} // namespace ttk
