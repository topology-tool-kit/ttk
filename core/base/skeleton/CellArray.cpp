#include <CellArray.h>

#include <iostream>

using namespace ttk;

CellArray::CellArray(const LongSimplexId *cellArray,
                     const LongSimplexId nbCells,
                     const unsigned char dimension)
  : cellArray_{cellArray}, nbCells_{nbCells}, dimension_{dimension} {
}

SimplexId CellArray::getCellVertexNumber(const LongSimplexId cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(cellId >= this->nbCells_) {
    std::cerr << "TTK: access to cell " << cellId << " on " << this->nbCells_
              << std::endl;
  }
#endif
  // Assume regular mesh
  return this->dimension_ + 1;
}

LongSimplexId CellArray::getCellVertex(const LongSimplexId cellId,
                                       const SimplexId localVertId) const {
  const SimplexId locNbVert = this->getCellVertexNumber(cellId);
#ifndef TTK_ENABLE_KAMIKAZE
  if(localVertId >= locNbVert) {
    std::cerr << "TTK: access to local vert " << localVertId << " on "
              << locNbVert << std::endl;
  }
#endif
  // Assume VTK < 9 layout
  return this->cellArray_[(locNbVert + 1) * cellId + 1 + localVertId];
}
