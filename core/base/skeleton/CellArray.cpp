#include <CellArray.h>

using namespace ttk;

#ifdef TTK_CELL_ARRAY_LEGACY
void CellArray::TranslateToFlatLayout(std::vector<LongSimplexId> &connectivity,
                                      std::vector<LongSimplexId> &offset,
                                      LongSimplexId *&singleArray) {
  const size_t coSize = connectivity.size();
  const size_t ofSize = offset.size();
  const size_t totalSize = coSize + ofSize - 1;
  singleArray = new LongSimplexId[totalSize];

  size_t con = 0;
  size_t single = 0;
  for(size_t off = 0; off < ofSize - 1; off++) {
    const LongSimplexId nbVerts = offset[off + 1] - offset[off];
    singleArray[single++] = nbVerts;
    for(int i = 0; i < nbVerts; i++) {
      singleArray[single++] = connectivity[con++];
    }
  }
}
#endif
