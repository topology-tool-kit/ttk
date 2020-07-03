#include <BarycentricSubdivision.h>

int ttk::BarycentricSubdivision::buildOutputTriangulation(
  ttk::ExplicitTriangulation &outputTriangl) {

  // ensure subdivision is already performed
  if(points_.empty() || cells_connectivity_.empty()) {
    return 1;
  }

  outputTriangl.setInputPoints(points_.size() / 3, points_.data());
#ifdef TTK_CELL_ARRAY_NEW
  outputTriangl.setInputCells(cells_offsets_.size() - 1,
                              cells_connectivity_.data(),
                              cells_offsets_.data());
#else
  // TODO use revert translator
  LongSimplexId *cells = nullptr;
  CellArray::TranslateToFlatLayout(cells_connectivity_, cells_offsets_, cells);
  outputTriangl.setInputCells(cells_offsets_.size() - 1, cells);
#endif

  return 0;
}
