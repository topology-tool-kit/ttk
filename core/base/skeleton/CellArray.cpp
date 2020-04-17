#include <CellArray.h>

using namespace ttk;

#ifdef CELL_ARRAY_NEW
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
