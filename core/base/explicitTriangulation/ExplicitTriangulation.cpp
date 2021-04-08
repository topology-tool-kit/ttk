#include <ExplicitTriangulation.h>
#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>
#include <cstring>

using namespace ttk;

ExplicitTriangulation::ExplicitTriangulation() {

  setDebugMsgPrefix("ExplicitTriangulation");

  clear();
}

ExplicitTriangulation::~ExplicitTriangulation() {
}

int ExplicitTriangulation::clear() {
  vertexNumber_ = 0;
  cellNumber_ = 0;
  doublePrecision_ = false;

  printMsg("Triangulation cleared.", debug::Priority::DETAIL);

  return AbstractTriangulation::clear();
}

size_t ExplicitTriangulation::footprint(size_t size) const {

  const auto printArrayFootprint
    = [this](const FlatJaggedArray &array, const std::string &name) {
        if(!array.empty() && !name.empty()) {
          this->printMsg(name + std::string{": "}
                         + std::to_string(array.footprint()) + " bytes");
        }
        return array.footprint();
      };

  size += printArrayFootprint(vertexNeighborData_, "vertexNeighborData_");
  size += printArrayFootprint(cellNeighborData_, "cellNeighborData_");
  size += printArrayFootprint(vertexEdgeData_, "vertexEdgeData_");
  size += printArrayFootprint(vertexTriangleData_, "vertexTriangleData_");
  size += printArrayFootprint(edgeTriangleData_, "edgeTriangleData_");
  size += printArrayFootprint(vertexStarData_, "vertexStarData_");
  size += printArrayFootprint(edgeStarData_, "edgeStarData_");
  size += printArrayFootprint(triangleStarData_, "triangleStarData_");
  size += printArrayFootprint(vertexLinkData_, "vertexLinkData_");
  size += printArrayFootprint(edgeLinkData_, "edgeLinkData_");
  size += printArrayFootprint(triangleLinkData_, "triangleLinkData_");

  return AbstractTriangulation::footprint(size);
}

int ExplicitTriangulation::preconditionBoundaryEdgesInternal() {
  if((!boundaryEdges_.empty()) && (boundaryEdges_.size() == edgeList_.size())) {
    return 0;
  }

  Timer tm{};

  preconditionEdgesInternal();
  boundaryEdges_.resize(edgeList_.size(), false);

  if(getDimensionality() == 2) {
    preconditionEdgeStarsInternal();
    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber(); i++) {
      if(edgeStarData_.size(i) == 1) {
        boundaryEdges_[i] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTriangleStarsInternal();
    preconditionTriangleEdgesInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        for(int j = 0; j < 3; j++) {
          boundaryEdges_[triangleEdgeList_[i][j]] = true;
        }
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary edges", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryTrianglesInternal() {
  if(getDimensionality() == 2)
    return 0;

  Timer tm{};

  if((!boundaryTriangles_.empty())
     && (boundaryTriangles_.size() == triangleList_.size())) {
    return 0;
  }

  preconditionTrianglesInternal();
  boundaryTriangles_.resize(triangleList_.size(), false);

  if(getDimensionality() == 3) {
    preconditionTriangleStarsInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        boundaryTriangles_[i] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary triangles", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryVerticesInternal() {
  if((!boundaryVertices_.empty())
     && ((SimplexId)boundaryVertices_.size() == vertexNumber_))
    return 0;

  Timer tm{};

  boundaryVertices_.resize(vertexNumber_, false);

  // create the list of boundary elements
  // create their star
  // look for singletons
  if(getDimensionality() == 1) {
    preconditionVertexStarsInternal();
    for(size_t i = 0; i < vertexStarData_.subvectorsNumber(); i++) {
      if(vertexStarData_.size(i) == 1) {
        boundaryVertices_[i] = true;
      }
    }
  } else if(getDimensionality() == 2) {
    preconditionEdgesInternal();
    preconditionEdgeStarsInternal();

    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber(); i++) {
      if(edgeStarData_.size(i) == 1) {
        boundaryVertices_[edgeList_[i][0]] = true;
        boundaryVertices_[edgeList_[i][1]] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTrianglesInternal();
    preconditionTriangleStarsInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        boundaryVertices_[triangleList_[i][0]] = true;
        boundaryVertices_[triangleList_[i][1]] = true;
        boundaryVertices_[triangleList_[i][2]] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary vertices", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionCellEdgesInternal() {
  OneSkeleton os;
  os.setWrapper(this);

  if(tetraEdgeList_.empty() && getDimensionality() == 3) {
    os.buildEdgeList(
      vertexNumber_, *cellArray_, nullptr, nullptr, &tetraEdgeList_);
  } else if(triangleEdgeList_.empty() && getDimensionality() == 2) {
    os.buildEdgeList(
      vertexNumber_, *cellArray_, nullptr, nullptr, &triangleEdgeList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellNeighborsInternal() {

  if(cellNeighborData_.empty()) {
    if(getDimensionality() == 3) {
      ThreeSkeleton threeSkeleton;
      threeSkeleton.setWrapper(this);
      threeSkeleton.buildCellNeighborsFromTriangles(
        vertexNumber_, *cellArray_, cellNeighborData_, &triangleStarData_);
    } else if(getDimensionality() == 2) {
      TwoSkeleton twoSkeleton;
      twoSkeleton.setWrapper(this);
      twoSkeleton.buildCellNeighborsFromEdges(
        vertexNumber_, *cellArray_, cellNeighborData_, &edgeStarData_);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellTrianglesInternal() {

  if(!tetraTriangleList_.size()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    if(triangleList_.size()) {
      // we already computed this guy, let's just get the cell triangles
      if(!triangleStarData_.empty()) {
        return twoSkeleton.buildTriangleList(
          vertexNumber_, *cellArray_, nullptr, nullptr, &tetraTriangleList_);
      } else {
        // let's compute the triangle star while we're at it...
        // it's just a tiny overhead.
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             nullptr, &triangleStarData_,
                                             &tetraTriangleList_);
      }
    } else {
      // we have not computed this guy, let's do it while we're at it
      if(!triangleStarData_.empty()) {
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             &triangleList_, nullptr,
                                             &tetraTriangleList_);
      } else {
        // let's compute the triangle star while we're at it...
        // it's just a tiny overhead.
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             &triangleList_, &triangleStarData_,
                                             &tetraTriangleList_);
      }
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgesInternal() {

  if(!edgeList_.size()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    // also computes edgeStar and triangleEdge / tetraEdge lists for free...
    if(getDimensionality() == 2) {
      return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_,
                                       &edgeStarData_, &triangleEdgeList_);
    } else if(getDimensionality() == 3) {
      return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_,
                                       &edgeStarData_, &tetraEdgeList_);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgeLinksInternal() {

  if(edgeLinkData_.empty()) {

    if(getDimensionality() == 2) {
      preconditionEdgesInternal();
      preconditionEdgeStarsInternal();

      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeLinks(
        edgeList_, edgeStarData_, *cellArray_, edgeLinkData_);
    } else if(getDimensionality() == 3) {
      preconditionEdgesInternal();
      preconditionEdgeStarsInternal();
      preconditionCellEdgesInternal();

      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeLinks(
        edgeList_, edgeStarData_, tetraEdgeList_, edgeLinkData_);
    } else {
      // unsupported dimension
      printErr("Unsupported dimension for edge link precondition");
      return -1;
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgeStarsInternal() {

  if(edgeStarData_.empty()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    return oneSkeleton.buildEdgeList<3>(
      vertexNumber_, *cellArray_, nullptr, &edgeStarData_, nullptr);
  }
  return 0;
}

int ExplicitTriangulation::preconditionEdgeTrianglesInternal() {

  if(edgeTriangleData_.empty()) {
    preconditionTriangleEdgesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildEdgeTriangles(vertexNumber_, *cellArray_,
                                          edgeTriangleData_, &edgeList_,
                                          &triangleEdgeList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTrianglesInternal() {

  if(!triangleList_.size()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_, &triangleList_,
                                  &triangleStarData_, &tetraTriangleList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleEdgesInternal() {

  if(!triangleEdgeList_.size()) {

    // WARNING
    // here triangleStarList and cellTriangleList will be computed (for
    // free) although they are not requireed to get the edgeTriangleList.
    // if memory usage is an issue, please change these pointers by nullptr.

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    return twoSkeleton.buildTriangleEdgeList(
      vertexNumber_, *cellArray_, triangleEdgeList_, &vertexEdgeData_,
      &edgeList_, &triangleList_, &triangleStarData_, &tetraTriangleList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleLinksInternal() {

  if(triangleLinkData_.empty()) {

    preconditionTriangleStarsInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildTriangleLinks(
      triangleList_, triangleStarData_, *cellArray_, triangleLinkData_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleStarsInternal() {

  if(triangleStarData_.empty()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildTriangleList(
      vertexNumber_, *cellArray_, &triangleList_, &triangleStarData_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionVertexEdgesInternal() {

  if((SimplexId)vertexEdgeData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;

    if(!edgeList_.size()) {
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_);
    }

    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexEdges(
      vertexNumber_, edgeList_, vertexEdgeData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexLinksInternal() {

  if((SimplexId)vertexLinkData_.subvectorsNumber() != vertexNumber_) {

    if(getDimensionality() == 2) {
      preconditionEdgesInternal();
      preconditionVertexStarsInternal();

      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(
        vertexStarData_, triangleEdgeList_, edgeList_, vertexLinkData_);
    } else if(getDimensionality() == 3) {
      preconditionTrianglesInternal();
      preconditionVertexStarsInternal();

      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(
        vertexStarData_, tetraTriangleList_, triangleList_, vertexLinkData_);
    } else {
      // unsupported dimension
      printErr("Unsupported dimension for vertex link precondition");
      return -1;
    }
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexNeighborsInternal() {

  if((SimplexId)vertexNeighborData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexNeighbors(
      vertexNumber_, *cellArray_, vertexNeighborData_, &edgeList_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexStarsInternal() {

  if((SimplexId)vertexStarData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);

    return zeroSkeleton.buildVertexStars(
      vertexNumber_, *cellArray_, vertexStarData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexTrianglesInternal() {

  if((SimplexId)vertexTriangleData_.subvectorsNumber() != vertexNumber_) {

    preconditionTrianglesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildVertexTriangles(
      vertexNumber_, triangleList_, vertexTriangleData_);
  }

  return 0;
}

template <typename T>
void writeBin(std::ofstream &stream, const T var) {
  stream.write(reinterpret_cast<const char *>(&var), sizeof(var));
}

// initialize static member variables
const char *ExplicitTriangulation::magicBytes_ = "TTKTriangulationFileFormat";
const unsigned long ExplicitTriangulation::formatVersion_ = 1;

int ExplicitTriangulation::writeToFile(std::ofstream &stream) const {

  // 1. magic bytes (char *)
  stream.write(this->magicBytes_, std::strlen(this->magicBytes_));
  // 2. format version (unsigned long)
  writeBin(stream, this->formatVersion_);
  // 3. dimensionality (int)
  const auto dim = this->getDimensionality();
  writeBin(stream, dim);
  // 4. number of vertices (SimplexId)
  const auto nVerts = this->getNumberOfVertices();
  writeBin(stream, nVerts);
  // 5. number of edges (SimplexId)
  const auto nEdges = this->getNumberOfEdges();
  writeBin(stream, nEdges);
  // 6. number of triangles (SimplexId, 0 in 1D)
  const auto nTriangles = dim > 1 ? this->getNumberOfTriangles() : 0;
  writeBin(stream, nTriangles);
  // 7. number of tetrahedron (SimplexId, 0 in 2D)
  const auto nTetras = dim > 2 ? this->getNumberOfCells() : 0;
  writeBin(stream, nTetras);

  // only write buffers oustside this->cellArray_ (cellVertex, vertexCoords),
  // those ones will be provided by VTK

  // fixed-size arrays (in AbstractTriangulation.h)

#define WRITE_FIXED(ARRAY)      \
  for(const auto &it : ARRAY) { \
    for(const auto el : it) {   \
      writeBin(stream, el);     \
    }                           \
  }

  // 8. edgeList (SimplexId array)
  WRITE_FIXED(this->edgeList_);
  // 9. triangleList (SimplexId array)
  WRITE_FIXED(this->triangleList_);
  // 10. triangleEdgeList (SimplexId array)
  WRITE_FIXED(this->triangleEdgeList_);
  // 11. tetraEdgeList (SimplexId array)
  WRITE_FIXED(this->tetraEdgeList_);
  // 12. tetraTriangleList (SimplexId array)
  WRITE_FIXED(this->tetraTriangleList_);

  // variable-size arrays (FlatJaggedArray in ExplicitTriangulation.h)

  const auto write_variable = [&stream](const FlatJaggedArray &arr) {
    for(size_t i = 0; i < arr.subvectorsNumber() + 1; ++i) {
      writeBin(stream, arr.offset(i));
    }
    for(size_t i = 0; i < arr.subvectorsNumber(); ++i) {
      const auto s = arr.size(i);
      for(SimplexId j = 0; j < s; ++j) {
        writeBin(stream, arr.get(i, j));
      }
    }
  };

  // 13. vertexNeighbors (SimplexId array, offsets then data)
  write_variable(this->vertexNeighborData_);
  // 14. cellNeighbors (SimplexId array, offsets then data)
  write_variable(this->cellNeighborData_);
  // 15. vertexEdges (SimplexId array, offsets then data)
  write_variable(this->vertexEdgeData_);
  // 16. edgeTriangles (SimplexId array, offsets then data)
  write_variable(this->edgeTriangleData_);
  // 17. vertexStars (SimplexId array, offsets then data)
  write_variable(this->vertexStarData_);
  // 18. edgeStars (SimplexId array, offsets then data)
  write_variable(this->edgeStarData_);
  // 19. triangleStars (SimplexId array, offsets then data)
  write_variable(this->triangleStarData_);
  // 20. vertexLinks (SimplexId array, offsets then data)
  write_variable(this->vertexLinkData_);
  // 21. edgeLinks (SimplexId array, offsets then data)
  write_variable(this->edgeLinkData_);
  // 22. triangleLinks (SimplexId array, offsets then data)
  write_variable(this->triangleLinkData_);

  // 23. boundary vertices (bool array)
  for(SimplexId i = 0; i < nVerts; ++i) {
    writeBin(stream, this->boundaryVertices_[i]);
  }
  // 24. boundary edges (bool array)
  for(SimplexId i = 0; i < nEdges; ++i) {
    writeBin(stream, this->boundaryEdges_[i]);
  }
  // 25. boundary triangles (bool array)
  for(SimplexId i = 0; i < nTriangles; ++i) {
    writeBin(stream, this->boundaryTriangles_[i]);
  }

  return 0;
}

template <typename T>
void readBin(std::ifstream &stream, T &res) {
  stream.read(reinterpret_cast<char *>(&res), sizeof(res));
}

template <typename T>
void readBinArray(std::ifstream &stream, T *res, size_t size) {
  stream.read(reinterpret_cast<char *>(res), size * sizeof(T));
}

int ExplicitTriangulation::readFromFile(std::ifstream &stream) {

  // 1. magic bytes (char *)
  const auto magicBytesLen = std::strlen(this->magicBytes_);
  std::vector<char> mBytes(magicBytesLen + 1);
  stream.read(mBytes.data(), magicBytesLen);
  const auto hasMagicBytes = std::strcmp(mBytes.data(), this->magicBytes_) == 0;
  if(!hasMagicBytes) {
    this->printErr("Could not find magic bytes in input files!");
    this->printErr("Aborting...");
    return 0;
  }
  // 2. format version (unsigned long)
  unsigned long version{};
  readBin(stream, version);
  if(version != this->formatVersion_) {
    this->printWrn("File format version (" + std::to_string(version)
                   + ") and software version ("
                   + std::to_string(this->formatVersion_) + ") are different!");
  }

  int dim{};
  SimplexId nVerts{}, nEdges{}, nTriangles{}, nTetras{};

  // 3. dimensionality (int)
  readBin(stream, dim);
  // 4. number of vertices (SimplexId)
  readBin(stream, nVerts);
  // 5. number of edges (SimplexId)
  readBin(stream, nEdges);
  // 6. number of triangles (SimplexId, 0 in 1D)
  readBin(stream, nTriangles);
  // 7. number of tetrahedron (SimplexId, 0 in 2D)
  readBin(stream, nTetras);

  if(dim != this->getDimensionality()) {
    this->printErr("Incorrect dimension!");
    return 0;
  }
  if(nVerts != this->getNumberOfVertices()) {
    this->printErr("Incorrect number of vertices!");
    return 0;
  }
  if((dim == 2 && nTriangles != this->getNumberOfCells())
     || (dim == 3 && nTetras != this->getNumberOfCells())) {
    this->printErr("Incorrect number of cells!");
    return 0;
  }

  // resize buffers
  this->edgeList_.resize(nEdges);
  this->triangleList_.resize(nTriangles);
  this->triangleEdgeList_.resize(nTriangles);
  this->tetraEdgeList_.resize(nTetras);
  this->tetraTriangleList_.resize(nTetras);

  // fixed-size arrays (in AbstractTriangulation.h)

  // 8. edgeList (SimplexId array)
  readBinArray(stream, this->edgeList_.data(), nEdges);
  // 9. triangleList (SimplexId array)
  readBinArray(stream, this->triangleList_.data(), nTriangles);
  // 10. triangleEdgeList (SimplexId array)
  readBinArray(stream, this->triangleEdgeList_.data(), nTriangles);
  // 11. tetraEdgeList (SimplexId array)
  readBinArray(stream, this->tetraEdgeList_.data(), nTetras);
  // 12. tetraTriangleList (SimplexId array)
  readBinArray(stream, this->tetraTriangleList_.data(), nTetras);

  // variable-size arrays (FlagJaggedArrays in ExplicitTriangulation.h)

  const auto read_variable
    = [&stream](FlatJaggedArray &arr, const SimplexId n_items) {
        std::vector<SimplexId> offsets{}, data{};
        offsets.resize(n_items + 1);
        readBinArray(stream, offsets.data(), offsets.size());
        data.resize(offsets.back());
        readBinArray(stream, data.data(), data.size());
        arr.setData(std::move(data), std::move(offsets));
      };

  // 13. vertexNeighbors (SimplexId array, offsets then data)
  read_variable(this->vertexNeighborData_, nVerts);
  // 14. cellNeighbors (SimplexId array, offsets then data)
  read_variable(this->cellNeighborData_, this->getNumberOfCells());
  // 15. vertexEdges (SimplexId array, offsets then data)
  read_variable(this->vertexEdgeData_, nVerts);
  // 16. edgeTriangles (SimplexId array, offsets then data)
  read_variable(this->edgeTriangleData_, nEdges);
  // 17. vertexStars (SimplexId array, offsets then data)
  read_variable(this->vertexStarData_, nVerts);
  // 18. edgeStars (SimplexId array, offsets then data)
  read_variable(this->edgeStarData_, nEdges);
  // 19. triangleStars (SimplexId array, offsets then data)
  read_variable(this->triangleStarData_, nTriangles);
  // 20. vertexLinks (SimplexId array, offsets then data)
  read_variable(this->vertexLinkData_, nVerts);
  // 21. edgeLinks (SimplexId array, offsets then data)
  read_variable(this->edgeLinkData_, nEdges);
  // 22. triangleLinks (SimplexId array, offsets then data)
  read_variable(this->triangleLinkData_, nTriangles);

  // resize more buffers
  this->boundaryVertices_.resize(nVerts);
  this->boundaryEdges_.resize(nEdges);
  this->boundaryTriangles_.resize(nTriangles);

  const auto read_bool
    = [&stream](std::vector<bool> &arr, const SimplexId n_items) {
        for(SimplexId i = 0; i < n_items; ++i) {
          char b{};
          stream.read(&b, sizeof(b));
          arr[i] = static_cast<bool>(b);
        }
      };

  // 23. boundary vertices (bool array)
  read_bool(this->boundaryVertices_, nVerts);
  // 24. boundary edges (bool array)
  read_bool(this->boundaryEdges_, nEdges);
  // 25. boundary triangles (bool array)
  read_bool(this->boundaryTriangles_, nTriangles);

  return 0;
}
