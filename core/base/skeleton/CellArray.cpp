#include <CellArray.h>

#include <iostream>

using namespace ttk;

CellArray::CellArray(const LongSimplexId *cellArray,
                     const LongSimplexId nbCells,
                     const unsigned char dimension)
  : cellArray_{cellArray}, nbCells_{nbCells}, dimension_{dimension} {
}
