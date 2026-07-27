/// \ingroup base
/// \class ttk::CellArray
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date March 2020.
///
/// \brief %CellArray generic array of cells
///
/// %CellArray is a generic container that allows to deal with various cell
/// layouts in memory \sa Triangulation \sa ttkTriangulation
/// This version assumes the cell data to be given by two
/// arrays:
/// * Connectivity, containing the list of point ids of each cells
/// * Offsets, containing the beginning and end position of each cell in the
/// Connectivity array See VTK 9 vtkCellArray documentation for more details
/// about this layout.
/// This class is more a view on data than a complete container. It does
/// not support memory modification operations.
/// This file will be installed by CMake as CellArray.h if the
/// TTK_CELL_ARRAY_LAYOUT is set to OffsetAndConnectivity

#ifndef _CELLARRAY_H
#define _CELLARRAY_H

#include <DataTypes.h>

#include <vector>

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

namespace ttk {
  // This is not a Debug class as we want to keep it as light as possible
  class CellArray {
  public:
    CellArray(const LongSimplexId *connectivity,
              const LongSimplexId *offset,
              const LongSimplexId nbCells)
      : connectivity_{connectivity}, offset_{offset}, nbCells_{nbCells} {
    }

    virtual ~CellArray() {
      if(ownerShip_) {
        delete[] connectivity_;
        delete[] offset_;
      }
    }

    /// Deal with data ownership
    void setOwnership(const bool o) {
      ownerShip_ = o;
    }

    /// Get the number of cells in the array
    inline LongSimplexId getNbCells() const {
      return nbCells_;
    }

    /// Get the number of vertices in the cell with the id: cellid
    /// Can deal with heterogeneous meshes
    /// \param cellId global id of the cell
    /// \return the offset difference between this cell and next.
    inline SimplexId getCellVertexNumber(const LongSimplexId cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(cellId >= this->nbCells_) {
        std::cerr << "TTK: access to cell " << cellId << " on "
                  << this->nbCells_ << std::endl;
      }
#endif
      // Connectivity size is nbCell + 1 for last cell
      return this->offset_[cellId + 1] - this->offset_[cellId];
    }
    /// Get the vertex id of the "localVertId"'nt vertex of the cell.
    /// Can deal with heterogeneous meshes
    /// \param cellId global id of the cell
    /// \param localVertId id of the vertex local to the cell
    /// \return the global id of the vertex
    inline LongSimplexId getCellVertex(const LongSimplexId cellId,
                                       const SimplexId localVertId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      const SimplexId locNbVert = this->getCellVertexNumber(cellId);
      if(localVertId >= locNbVert) {
        std::cerr << "TTK: access to local vert " << localVertId << " on "
                  << locNbVert << std::endl;
      }
#endif
      return connectivity_[offset_[cellId] + localVertId];
    }

  protected:
    const LongSimplexId *connectivity_;
    const LongSimplexId *offset_;
    const LongSimplexId nbCells_;
    bool ownerShip_ = false;
  };
} // namespace ttk

#endif
