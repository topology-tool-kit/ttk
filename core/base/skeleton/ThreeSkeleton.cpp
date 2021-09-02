#include <ThreeSkeleton.h>
#include <boost/container/small_vector.hpp>

using namespace ttk;

ThreeSkeleton::ThreeSkeleton() {
  // TODO could use typeid for automatic value
  setDebugMsgPrefix("ThreeSkeleton");
}

int ThreeSkeleton::buildCellNeighborsFromTriangles(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  FlatJaggedArray &cellNeighbors,
  FlatJaggedArray *triangleStars) const {

  auto localTriangleStars = triangleStars;
  FlatJaggedArray defaultTriangleStars{};
  if(!localTriangleStars) {
    localTriangleStars = &defaultTriangleStars;
  }

  if(localTriangleStars->empty()) {
    TwoSkeleton twoSkeleton;
    twoSkeleton.setThreadNumber(threadNumber_);
    twoSkeleton.setDebugLevel(debugLevel_);
    twoSkeleton.buildTriangleList(
      vertexNumber, cellArray, nullptr, localTriangleStars);
  }

  Timer t;

  printMsg("Building cell neighbors", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  const SimplexId triangleNumber = localTriangleStars->subvectorsNumber();
  std::vector<SimplexId> offsets(cellNumber + 1);
  // number of neighbors processed per cell
  std::vector<SimplexId> neighborsId(cellNumber);

  for(SimplexId i = 0; i < triangleNumber; i++) {
    if(localTriangleStars->size(i) == 2) {
      // tetra cells in triangle i's star
      const auto cs0 = localTriangleStars->get(i, 0);
      const auto cs1 = localTriangleStars->get(i, 1);
      offsets[cs0 + 1]++;
      offsets[cs1 + 1]++;
    }
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> neighbors(offsets.back());

  // fill flat neighbors vector using offsets and neighbors count vectors
  for(SimplexId i = 0; i < triangleNumber; i++) {
    if(localTriangleStars->size(i) == 2) {
      // tetra cells in triangle i's star
      const auto cs0 = localTriangleStars->get(i, 0);
      const auto cs1 = localTriangleStars->get(i, 1);
      neighbors[offsets[cs0] + neighborsId[cs0]] = cs1;
      neighborsId[cs0]++;
      neighbors[offsets[cs1] + neighborsId[cs1]] = cs0;
      neighborsId[cs1]++;
    }
  }

  // fill FlatJaggedArray struct
  cellNeighbors.setData(std::move(neighbors), std::move(offsets));

  printMsg("Built " + std::to_string(cellNumber) + " cell neighbors", 1,
           t.getElapsedTime(), 1);

  // ethaneDiol.vtu, 8.7Mtets, vger (4coresHT)
  // 1 thread: 9.80 s
  // 4 threads: 14.18157 s

  // ethaneDiol.vtu, 8.7Mtets, richard (8coresHT)
  // 1 thread: 9.6763 s
  // 4 threads: 16.7489 s

  // ethaneDiol.vtu, 8.7Mtets, vger (4coresHT), VTK implementation
  // 1 thread: 14.9023 s
  // 4 threads: 8.58048 s

  // ethaneDiol.vtu, 8.7Mtets, richard (8coresHT), VTK implementation
  // 1 thread: 15.1387  s
  // 8 threads: 5.38286 s

  return 0;
}

int ThreeSkeleton::buildCellNeighborsFromVertices(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  FlatJaggedArray &cellNeighbors,
  FlatJaggedArray *vertexStars) const {

  // TODO: ASSUME uniform mesh here!
  if(cellArray.getNbCells() && cellArray.getCellVertexNumber(0) == 3) {
    TwoSkeleton twoSkeleton;
    twoSkeleton.setDebugLevel(debugLevel_);
    twoSkeleton.setThreadNumber(threadNumber_);
    return twoSkeleton.buildCellNeighborsFromVertices(
      vertexNumber, cellArray, cellNeighbors, vertexStars);
  }

  if(cellArray.getNbCells() && cellArray.getCellVertexNumber(0) <= 2) {
    // 1D
    printErr("buildCellNeighborsFromVertices in 1D:");
    printErr("Not implemented! TODO?!");
    return -1;
  }

  auto localVertexStars = vertexStars;
  FlatJaggedArray defaultVertexStars{};

  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if(localVertexStars->empty()) {

    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.buildVertexStars(vertexNumber, cellArray, *localVertexStars);
  }

  Timer t;

  printMsg("Building cell neighnors", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  using boost::container::small_vector;
  // for each cell/tetra, a vector of neighbors
  std::vector<small_vector<SimplexId, 4>> neighbors(cellNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);

    // go triangle by triangle
    for(SimplexId j = 0; j < nbVertCell; j++) {

      SimplexId v0 = cellArray.getCellVertex(cid, j);
      SimplexId v1 = cellArray.getCellVertex(cid, (j + 1) % nbVertCell);
      SimplexId v2 = cellArray.getCellVertex(cid, (j + 2) % nbVertCell);

      // perform an intersection of the 3 (sorted) star lists
      SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
      SimplexId intersection = -1;

      while(pos0 < localVertexStars->size(v0)
            && pos1 < localVertexStars->size(v1)
            && pos2 < localVertexStars->size(v2)) {

        SimplexId biggest = localVertexStars->get(v0, pos0);
        if(localVertexStars->get(v1, pos1) > biggest) {
          biggest = localVertexStars->get(v1, pos1);
        }
        if(localVertexStars->get(v2, pos2) > biggest) {
          biggest = localVertexStars->get(v2, pos2);
        }

        for(SimplexId l = pos0; l < localVertexStars->size(v0); l++) {
          if(localVertexStars->get(v0, l) < biggest) {
            pos0++;
          } else {
            break;
          }
        }
        for(SimplexId l = pos1; l < localVertexStars->size(v1); l++) {
          if(localVertexStars->get(v1, l) < biggest) {
            pos1++;
          } else {
            break;
          }
        }
        for(SimplexId l = pos2; l < localVertexStars->size(v2); l++) {
          if(localVertexStars->get(v2, l) < biggest) {
            pos2++;
          } else {
            break;
          }
        }

        if(pos0 < localVertexStars->size(v0)
           && pos1 < localVertexStars->size(v1)
           && pos2 < localVertexStars->size(v2)) {

          if((localVertexStars->get(v0, pos0)
              == localVertexStars->get(v1, pos1))
             && (localVertexStars->get(v0, pos0)
                 == localVertexStars->get(v2, pos2))) {

            if(localVertexStars->get(v0, pos0) != cid) {
              intersection = localVertexStars->get(v0, pos0);
              break;
            }

            pos0++;
            pos1++;
            pos2++;
          }
        }
      }

      if(intersection != -1) {
        neighbors[cid].emplace_back(intersection);
      }
    }
  }

  // convert to a FlatJaggedArray
  cellNeighbors.fillFrom(neighbors);

  printMsg("Built " + std::to_string(cellNumber) + " cell neighbors", 1,
           t.getElapsedTime(), threadNumber_);

  // ethaneDiol.vtu, 8.7Mtets, richard (4coresHT)
  // 1 thread: 9.39488 s
  // 2 threads: 5.93132 s [faster than any other implementation]
  // 4 threads: 3.3445 s [~x3, not too bad]

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 4.4475 s [speedup on VTK ~2.5]
  // 24 threads: 1.4488 s

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 29.4951 s [speedup on VTK ~3.3]
  // 24 threads: 12.7058 s (parallel speed up ~3, speed up on VTK: 2)

  // NOTE:
  // when dealing with simple types such as int, it's better to have heap
  // variable (instead of a vector<int> indexed by threadId).
  // it is MUCH MUCH MUCH better
  // especially with write access (cache misses)

  return 0;
}
