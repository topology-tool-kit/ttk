#include <Geometry.h>
#include <Laplacian.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>

#include <array>

template <typename T, typename SparseMatrixType = Eigen::SparseMatrix<T>>
int ttk::Laplacian::discreteLaplacian(SparseMatrixType &output,
                                      const Triangulation &triangulation) {

  using Triplet = Eigen::Triplet<T>;
  auto vertexNumber = triangulation.getNumberOfVertices();
  auto edgeNumber = triangulation.getNumberOfEdges();

  // early return when input graph is empty
  if(vertexNumber <= 0) {
    return -1;
  }

#ifdef TTK_ENABLE_OPENMP
  // get thread number from triangulation?
  auto threadNumber = triangulation.getThreadNumber();
#endif // TTK_ENABLE_OPENMP

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
    SimplexId nneigh = triangulation.getVertexNeighborNumber(SimplexId(i));
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

  output.setFromTriplets(triplets.begin(), triplets.end());

  return 0;
}

template <typename T, typename SparseMatrixType = Eigen::SparseMatrix<T>>
int ttk::Laplacian::cotanWeights(SparseMatrixType &output,
                                 const Triangulation &triangulation) {

  using Triplet = Eigen::Triplet<T>;
  auto vertexNumber = triangulation.getNumberOfVertices();
  auto edgeNumber = triangulation.getNumberOfEdges();

#ifdef TTK_ENABLE_OPENMP
  // get thread number from triangulation?
  auto threadNumber = triangulation.getThreadNumber();
#endif // TTK_ENABLE_OPENMP

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

  // hold sum of cotan weights for every vertex
  std::vector<T> vertexWeightSum(vertexNumber, T(0.0));

  // iterate over all edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber; ++i) {

    // the two vertices of the current edge (+ a third)
    std::vector<SimplexId> edgeVertices(3);
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation.getEdgeVertex(i, j, edgeVertices[j]);
    }

    // get the triangles that share the current edge
    // in 2D only 2, in 3D, maybe more...
    SimplexId trianglesNumber = triangulation.getEdgeTriangleNumber(i);
    // stores the triangles ID for every triangle around the current edge
    std::vector<SimplexId> edgeTriangles(trianglesNumber);
    for(SimplexId j = 0; j < trianglesNumber; ++j) {
      triangulation.getEdgeTriangle(i, j, edgeTriangles[j]);
    }

    // iterate over current edge triangles
    std::vector<T> angles{};

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

    // store the cotan weight sum for the two vertices of the current edge
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
    {
      vertexWeightSum[edgeVertices[0]] += cotan_weight;
      vertexWeightSum[edgeVertices[1]] += cotan_weight;
    }
  }

  // on the diagonal: sum of cotan weights for every vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    triplets[2 * edgeNumber + i] = Triplet(i, i, vertexWeightSum[i]);
  }

  output.setFromTriplets(triplets.begin(), triplets.end());

  return 0;
}

// explicit intantiations for floating-point types
template int
  ttk::Laplacian::discreteLaplacian<float>(Eigen::SparseMatrix<float> &output,
                                           const Triangulation &triangulation);
template int
  ttk::Laplacian::discreteLaplacian<double>(Eigen::SparseMatrix<double> &output,
                                            const Triangulation &triangulation);

template int
  ttk::Laplacian::cotanWeights<float>(Eigen::SparseMatrix<float> &output,
                                      const Triangulation &triangulation);
template int
  ttk::Laplacian::cotanWeights<double>(Eigen::SparseMatrix<double> &output,
                                       const Triangulation &triangulation);

#endif // TTK_ENABLE_EIGEN
