#include <Geometry.h>
#include <Laplacian.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>

#include <array>

template <typename T,
          class TriangulationType,
          typename SparseMatrixType = Eigen::SparseMatrix<T>>
int ttk::Laplacian::discreteLaplacian(SparseMatrixType &output,
                                      const Debug &dbg,
                                      const TriangulationType &triangulation) {

  Timer tm{};

  using Triplet = Eigen::Triplet<T>;
  const auto vertexNumber = triangulation.getNumberOfVertices();
  const auto edgeNumber = triangulation.getNumberOfEdges();

  // early return when input graph is empty
  if(vertexNumber <= 0) {
    return -1;
  }

  const auto threadNumber = dbg.getThreadNumber();

  // clear output
  output.resize(vertexNumber, vertexNumber);
  output.setZero();

  // number of triplets to insert into laplacian matrix: vertexNumber_
  // values on the diagonal + 2 values per edge
  std::vector<Triplet> triplets(vertexNumber + 2 * edgeNumber);

  // on the diagonal: number of neighbors
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    const auto nneigh = triangulation.getVertexNeighborNumber(i);
    triplets[i] = Triplet(i, i, T(nneigh));
  }

  // neighbors mapping: loop over edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber; ++i) {
    // the two vertices of the current edge
    std::vector<SimplexId> edgeVertices(2);
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation.getEdgeVertex(i, j, edgeVertices[j]);
    }
    // fill triplets for both vertices of current edge
    triplets[vertexNumber + 2 * i]
      = Triplet(edgeVertices[0], edgeVertices[1], T(-1.0));
    triplets[vertexNumber + 2 * i + 1]
      = Triplet(edgeVertices[1], edgeVertices[0], T(-1.0));
  }

#ifndef __clang_analyzer__
  output.setFromTriplets(triplets.begin(), triplets.end());
#endif // __clang_analyzer__

  dbg.printMsg(
    "Computed Discrete Laplacian", 1.0, tm.getElapsedTime(), threadNumber);

  return 0;
}

template <typename T,
          class TriangulationType,
          typename SparseMatrixType = Eigen::SparseMatrix<T>>
int ttk::Laplacian::cotanWeights(SparseMatrixType &output,
                                 const Debug &dbg,
                                 const TriangulationType &triangulation) {

  Timer tm{};

  using Triplet = Eigen::Triplet<T>;
  const auto vertexNumber = triangulation.getNumberOfVertices();
  const auto edgeNumber = triangulation.getNumberOfEdges();

  const auto threadNumber = dbg.getThreadNumber();

  // early return when input graph is empty
  if(vertexNumber <= 0) {
    return -1;
  }

  // clear output
  output.resize(vertexNumber, vertexNumber);
  output.setZero();

  // number of triplets to insert into laplacian matrix: vertexNumber_
  // values on the diagonal + 2 values per edge
  std::vector<Triplet> triplets(vertexNumber + 2 * edgeNumber);

  std::vector<SimplexId> edgeTriangles{};
  std::vector<T> angles{};

  // iterate over all edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber) \
  firstprivate(edgeTriangles, angles)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber; ++i) {

    edgeTriangles.clear();
    angles.clear();

    // the two vertices of the current edge (+ a third)
    std::array<SimplexId, 3> edgeVertices{};
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation.getEdgeVertex(i, j, edgeVertices[j]);
    }

    // get the triangles that share the current edge
    // in 2D only 2, in 3D, maybe more...
    const auto trianglesNumber = triangulation.getEdgeTriangleNumber(i);
    // stores the triangles ID for every triangle around the current edge
    edgeTriangles.resize(trianglesNumber);
    for(SimplexId j = 0; j < trianglesNumber; ++j) {
      triangulation.getEdgeTriangle(i, j, edgeTriangles[j]);
    }

    // iterate over current edge triangles
    angles.reserve(trianglesNumber);

    for(const auto &j : edgeTriangles) {

      // get the third vertex of the triangle
      SimplexId thirdNeigh;
      // a triangle has only three vertices
      for(SimplexId k = 0; k < 3; ++k) {
        triangulation.getTriangleVertex(j, k, thirdNeigh);
        if(thirdNeigh != edgeVertices[0] && thirdNeigh != edgeVertices[1]) {
          // store the third vertex ID into the edgeVertices array to
          // be more easily handled
          edgeVertices[2] = thirdNeigh;
          break;
        }
      }
      // compute the 3D coords of the three vertices
      std::array<float, 9> coords{};
      for(SimplexId k = 0; k < 3; ++k) {
        triangulation.getVertexPoint(
          edgeVertices[k], coords[3 * k], coords[3 * k + 1], coords[3 * k + 2]);
      }
      angles.emplace_back(ttk::Geometry::angle(&coords[6], // edgeVertices[2]
                                               &coords[0], // edgeVertices[0]
                                               &coords[6], // edgeVertices[2]
                                               &coords[3]) // edgeVertices[1]
      );
    }

    // cotan weights for every triangle around the current edge
    T cotan_weight{0.0};
    // C++ has no map statement until C++17 (std::transform)
    for(auto &angle : angles) {
      cotan_weight += T(1.0) / std::tan(angle);
    }

    // since we iterate over the edges, fill the laplacian matrix
    // symmetrically for the two vertices
    triplets[2 * i] = Triplet(edgeVertices[0], edgeVertices[1], -cotan_weight);
    triplets[2 * i + 1]
      = Triplet(edgeVertices[1], edgeVertices[0], -cotan_weight);
  }

  // on the diagonal: sum of cotan weights for every vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    T vertWeightSum{};
    const auto nEdges{triangulation.getVertexEdgeNumber(i)};
    // get the (-) cotan weights from the edges triplets
    for(SimplexId j = 0; j < nEdges; ++j) {
      SimplexId e{};
      triangulation.getVertexEdge(i, j, e);
      vertWeightSum += triplets[2 * e].value();
    }

    triplets[2 * edgeNumber + i] = Triplet(i, i, -vertWeightSum);
  }

#ifndef __clang_analyzer__
  output.setFromTriplets(triplets.begin(), triplets.end());
#endif // __clang_analyzer__

  dbg.printMsg("Computed Laplacian with Cotan Weights", 1.0,
               tm.getElapsedTime(), threadNumber);

  return 0;
}

#define LAPLACIAN_SPECIALIZE(TYPE)                                         \
  template int ttk::Laplacian::discreteLaplacian<TYPE>(                    \
    Eigen::SparseMatrix<TYPE> &, const Debug &dbg, const Triangulation &); \
  template int ttk::Laplacian::cotanWeights<TYPE>(                         \
    Eigen::SparseMatrix<TYPE> &, const Debug &dbg, const Triangulation &)

// explicit intantiations for floating-point types
LAPLACIAN_SPECIALIZE(float);
LAPLACIAN_SPECIALIZE(double);

#endif // TTK_ENABLE_EIGEN
