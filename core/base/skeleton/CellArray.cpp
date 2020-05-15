#include <CellArray.h>

using namespace ttk;

#ifdef CELL_ARRAY_NEW
// TODO: Remove as it should not be used anymore
void CellArray::SingleToOffsetAndCo(const LongSimplexId *singleArray,
                                    size_t nbCells,
                                    std::vector<LongSimplexId> &connectivity,
                                    std::vector<LongSimplexId> &offset) {
  connectivity.reserve(nbCells * 3); // real size not known
  offset.resize(nbCells + 1);

  size_t curPos = 0;
  for(size_t cid = 0; cid < nbCells; cid++) {
    offset[cid] = connectivity.size();

    const LongSimplexId nbVerts = singleArray[curPos];
    curPos++;
    for(LongSimplexId v = 0; v < nbVerts; v++) {
      connectivity.emplace_back(singleArray[curPos]);
      curPos++;
    }
  }
  // last one
  offset[nbCells] = connectivity.size();
}
#endif

#ifdef CELL_ARRAY_LEGACY
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
