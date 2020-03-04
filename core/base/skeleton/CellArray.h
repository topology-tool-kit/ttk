/// \ingroup base
/// \class ttk::CellArray
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date March 2020.
///
/// \brief %CellArray generic array of cells
///
/// %CellArray is a generic container that allows to deal with various cell
/// layouts in memory \sa Triangulation \sa ttkTriangulation

#ifndef _CELLARRAY_H
#define _CELLARRAY_H

#include <DataTypes.h>

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

namespace ttk {

  // This is not a Debug class as we want to keep it as light as possible
  class CellArray {
  public:
    CellArray(const LongSimplexId *cellArray,
              const LongSimplexId nbCells,
              const unsigned char dimension);

    virtual ~CellArray() {
      if(ownerShip_) {
        delete[] cellArray_;
      }
    }

    /// Deal with data ownership
    void setOwnership(const bool o) {
      ownerShip_ = o;
    }

    /// Retrieve the dimension
    const unsigned char getDimension() const {
      return dimension_;
    }

    /// Get the number of cells in the array
    LongSimplexId getNbCells() const {
      return nbCells_;
    }

    /// Get the number of vertices in the cell with the id: cellid
    /// \param cellid global id of the cell
    /// \return dimension + 1 for now as we only accept regular meshes
    SimplexId getCellVertexNumber(const LongSimplexId cellid) const;

    /// Get the vertex id of the "localVertId"'nt vertex of the cell.
    /// \param cellId global id of the cell
    /// \param localVertId id of the vertex local to the cell,
    /// usually lower than 4 in 3D and lower than 3 in 2D.
    /// \return the global id of the vertex
    LongSimplexId getCellVertex(const LongSimplexId cellId,
                                const SimplexId localVertId) const;

  protected:
    const LongSimplexId *cellArray_;
    const LongSimplexId nbCells_;
    const unsigned char dimension_;
    bool ownerShip_ = false;
  };
} // namespace ttk

#endif
