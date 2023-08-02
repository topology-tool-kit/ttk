/// \ingroup base
/// \class ttk::CellArray
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date March 2020.
///
/// \brief %CellArray generic array of cells
///
/// %CellArray is a generic container that allows to deal with various cell
/// layouts in memory \sa Triangulation \sa ttkTriangulation
/// This version assumes the cell data to be given as a unique array.
/// This array contains for each cell the number of points in the cell then the
/// point ids of the cell. See VTK 8.2 vtkCellArray documentation for more
/// details about this layout.
/// This file will be installed by CMake as CellArray.h if the
/// TTK_CELL_ARRAY_LAYOUT is set to SingleArray

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
    CellArray(const LongSimplexId *cellArray,
              const LongSimplexId nbCells,
              const unsigned char dimension)
      : cellArray_{cellArray}, nbCells_{nbCells}, dimension_{dimension} {
    }

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
    inline unsigned char getDimension() const {
      return dimension_;
    }

    /// Get the number of cells in the array
    inline LongSimplexId getNbCells() const {
      return nbCells_;
    }

    /// Get the number of vertices in the cell with the id: cellid
    /// \param cellId global id of the cell
    /// \return dimension + 1 for now as we only accept regular meshes
    inline SimplexId getCellVertexNumber(const LongSimplexId cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(cellId >= this->nbCells_) {
        std::cerr << "TTK: access to cell " << cellId << " on "
                  << this->nbCells_ << std::endl;
      }
#endif
      // WARNING: ASSUME Regular Mesh here
      return this->dimension_ + 1;
    }
    /// Get the vertex id of the "localVertId"'nt vertex of the cell.
    /// \param cellId global id of the cell
    /// \param localVertId id of the vertex local to the cell,
    /// usually lower than 4 in 3D and lower than 3 in 2D.
    /// \return the global id of the vertex
    inline LongSimplexId getCellVertex(const LongSimplexId cellId,
                                       const SimplexId localVertId) const {
      const SimplexId locNbVert = this->getCellVertexNumber(cellId);
#ifndef TTK_ENABLE_KAMIKAZE
      if(localVertId >= locNbVert) {
        std::cerr << "TTK: access to local vert " << localVertId << " on "
                  << locNbVert << std::endl;
      }
#endif
      // Assume VTK < 9 layout and uniform mesh (only one type of cells)
      return this->cellArray_[(locNbVert + 1) * cellId + 1 + localVertId];
    }

    static void TranslateToFlatLayout(std::vector<LongSimplexId> &connectivity,
                                      std::vector<LongSimplexId> &offset,
                                      LongSimplexId *&singleArray);

  protected:
    const LongSimplexId *cellArray_;
    const LongSimplexId nbCells_;
    const unsigned char dimension_;
    bool ownerShip_ = false;
  };
} // namespace ttk

#endif
