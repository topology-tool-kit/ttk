/// \ingroup base
/// \class ttk::Icosphere
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// This filter creates an Icosphere with a specified radius, center, and number
/// of subdivisions.

#pragma once

// base code includes
#include <Debug.h>

// std includes
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <vector>

namespace ttk {

  class Icosphere : virtual public Debug {

  public:
    Icosphere() {
      this->setDebugMsgPrefix("Icosphere");
    }
    ~Icosphere() override = default;

    /**
     * Efficiently computes for a given subdivision level the number of
     * resulting vertices and triangles.
     */
    int computeNumberOfVerticesAndTriangles(size_t &nVertices,
                                            size_t &nTriangles,
                                            const size_t nSubdivisions) const {
      nVertices = 12;
      nTriangles = 20;

      for(size_t i = 0; i < nSubdivisions; i++) {
        nVertices += (nTriangles * 3) / 2;
        nTriangles = nTriangles * 4;
      }

      return 1;
    }

    /**
     * Translates an icosphere
     */
    template <typename DT, typename IT>
    int translateIcosphere(
      // Output
      DT *vertexCoords,
      IT *connectivityList,

      // Input
      const size_t &icosphereIndex,
      const size_t &nVerticesPerIcosphere,
      const size_t &nTrianglesPerIcosphere,
      const DT *centers) const;

    /**
     * Computes an icosphere for a given subdivision level.
     */
    template <typename DT, typename IT>
    int computeIcosphere(
      // Output
      DT *vertexCoords,
      IT *connectivityList,

      // Input
      const size_t &nSubdivisions) const;

    /**
     * Computes an icosphere for a given subdivision level, radius, and center.
     */
    template <typename DT, typename IT>
    int computeIcospheres(
      // Output
      DT *vertexCoords,
      IT *connectivityList,

      // Input
      const size_t &nSpheres,
      const size_t &nSubdivisions,
      const DT &radius,
      const DT *center,

      // Optional Output
      DT *normals = nullptr) const;

  private:
    /**
     * Adds the coordinates of a vertex to the vertexCoords array at
     * vertexIndex*3, return the current vertexIndex, and increases the
     * vertexIndex by one.
     */
    template <typename DT, typename IT>
    IT addVertex(const DT &x,
                 const DT &y,
                 const DT &z,
                 DT *vertexCoords,
                 IT &vertexIndex) const {
      IT cIndex = vertexIndex * 3;

      DT length = std::sqrt(x * x + y * y + z * z);
      vertexCoords[cIndex] = x / length;
      vertexCoords[cIndex + 1] = y / length;
      vertexCoords[cIndex + 2] = z / length;

      return vertexIndex++;
    }

    /**
     * Adds the ijk ids of a new triangle to the connectivityList at
     * triangleIndex*3, returns the current triangleIndex, and increases
     * the triangleIndex by one.
     */
    template <class IT>
    IT addTriangle(const IT &i,
                   const IT &j,
                   const IT &k,
                   IT *connectivityList,
                   IT &triangleIndex) const {
      IT cIndex = triangleIndex * 3;
      connectivityList[cIndex + 0] = i;
      connectivityList[cIndex + 1] = j;
      connectivityList[cIndex + 2] = k;
      return triangleIndex++;
    }

    /**
     * Adds a new vertex at the middle of an edge between the vertices i
     * and j. If this vertex was already created in a previous call (as
     * an edge of an already processed triangle), the function only
     * returns the index of the existing vertex.
     */
    template <typename DT, typename IT>
    IT addMidVertex(
      const IT &i,
      const IT &j,
      std::unordered_map<std::pair<IT, IT>, IT, boost::hash<std::pair<IT, IT>>>
        &processedEdges,
      DT *vertexCoords,
      IT &vertexIndex) const {
      bool firstIsSmaller = i < j;
      IT a = firstIsSmaller ? i : j;
      IT b = firstIsSmaller ? j : i;

      // Check if edge was already processed
      {
        auto it = processedEdges.find({a, b});
        if(it != processedEdges.end())
          return it->second;
      }

      // Otherwise add mid vertex
      {
        IT aOffset = a * 3;
        IT bOffset = b * 3;
        DT mx = (vertexCoords[aOffset] + vertexCoords[bOffset]) / 2.0;
        DT my = (vertexCoords[aOffset + 1] + vertexCoords[bOffset + 1]) / 2.0;
        DT mz = (vertexCoords[aOffset + 2] + vertexCoords[bOffset + 2]) / 2.0;

        IT mIndex
          = this->addVertex<DT, IT>(mx, my, mz, vertexCoords, vertexIndex);
        processedEdges.insert({{a, b}, mIndex});
        return mIndex;
      }
    }
  };
} // namespace ttk

template <typename DT, typename IT>
int ttk::Icosphere::computeIcosphere(
  // Output
  DT *vertexCoords,
  IT *connectivityList,

  // Input
  const size_t &nSubdivisions) const {

  // initialize timer and memory
  Timer timer;

  // print status
  const std::string msg
    = "Computing Icosphere (S: " + std::to_string(nSubdivisions) + ")";
  this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  IT vertexIndex = 0;
  IT triangleIndex = 0;

  // build icosahedron
  {
    // create 12 vertices
    DT t = (1.0 + sqrt(5.0)) / 2.0;
    this->addVertex<DT, IT>(-1, t, 0, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(1, t, 0, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(-1, -t, 0, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(1, -t, 0, vertexCoords, vertexIndex);

    this->addVertex<DT, IT>(0, -1, t, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(0, 1, t, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(0, -1, -t, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(0, 1, -t, vertexCoords, vertexIndex);

    this->addVertex<DT, IT>(t, 0, -1, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(t, 0, 1, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(-t, 0, -1, vertexCoords, vertexIndex);
    this->addVertex<DT, IT>(-t, 0, 1, vertexCoords, vertexIndex);

    // create 20 triangles
    this->addTriangle<IT>(0, 11, 5, connectivityList, triangleIndex);
    this->addTriangle<IT>(0, 5, 1, connectivityList, triangleIndex);
    this->addTriangle<IT>(0, 1, 7, connectivityList, triangleIndex);
    this->addTriangle<IT>(0, 7, 10, connectivityList, triangleIndex);
    this->addTriangle<IT>(0, 10, 11, connectivityList, triangleIndex);
    this->addTriangle<IT>(1, 5, 9, connectivityList, triangleIndex);
    this->addTriangle<IT>(5, 11, 4, connectivityList, triangleIndex);
    this->addTriangle<IT>(11, 10, 2, connectivityList, triangleIndex);
    this->addTriangle<IT>(10, 7, 6, connectivityList, triangleIndex);
    this->addTriangle<IT>(7, 1, 8, connectivityList, triangleIndex);
    this->addTriangle<IT>(3, 9, 4, connectivityList, triangleIndex);
    this->addTriangle<IT>(3, 4, 2, connectivityList, triangleIndex);
    this->addTriangle<IT>(3, 2, 6, connectivityList, triangleIndex);
    this->addTriangle<IT>(3, 6, 8, connectivityList, triangleIndex);
    this->addTriangle<IT>(3, 8, 9, connectivityList, triangleIndex);
    this->addTriangle<IT>(4, 9, 5, connectivityList, triangleIndex);
    this->addTriangle<IT>(2, 4, 11, connectivityList, triangleIndex);
    this->addTriangle<IT>(6, 2, 10, connectivityList, triangleIndex);
    this->addTriangle<IT>(8, 6, 7, connectivityList, triangleIndex);
    this->addTriangle<IT>(9, 8, 1, connectivityList, triangleIndex);
  }

  // refine icosahedron
  if(nSubdivisions > 0) {
    // create temporary connectivityList
    size_t nVertices = 0;
    size_t nTriangles = 0;
    this->computeNumberOfVerticesAndTriangles(
      nVertices, nTriangles, nSubdivisions);

    std::vector<IT> connectivityListTemp(nTriangles * 3, 0);

    // cache to store processed edges
    std::unordered_map<std::pair<IT, IT>, IT, boost::hash<std::pair<IT, IT>>>
      processedEdges;

    // iterate over nSubdivisions
    for(size_t s = 0; s < nSubdivisions; s++) {
      // swap lists
      const IT *oldList
        = s % 2 == 0 ? connectivityList : connectivityListTemp.data();
      IT *newList = s % 2 != 0 ? connectivityList : connectivityListTemp.data();

      // reset indicies
      const size_t nOldTriangles = triangleIndex;
      triangleIndex = 0;

      for(size_t i = 0; i < nOldTriangles; i++) {
        const IT offset = i * 3;

        // compute mid points
        const IT a
          = this->addMidVertex(oldList[offset + 0], oldList[offset + 1],
                               processedEdges, vertexCoords, vertexIndex);
        const IT b
          = this->addMidVertex(oldList[offset + 1], oldList[offset + 2],
                               processedEdges, vertexCoords, vertexIndex);
        const IT c
          = this->addMidVertex(oldList[offset + 2], oldList[offset + 0],
                               processedEdges, vertexCoords, vertexIndex);

        // replace triangle by 4 triangles
        this->addTriangle(oldList[offset + 0], a, c, newList, triangleIndex);
        this->addTriangle(oldList[offset + 1], b, a, newList, triangleIndex);
        this->addTriangle(oldList[offset + 2], c, b, newList, triangleIndex);
        this->addTriangle(a, b, c, newList, triangleIndex);
      }

      // print progress
      this->printMsg(
        msg, ((DT)s) / nSubdivisions, ttk::debug::LineMode::REPLACE);
    }

    // if uneven number of nSubdivisions then copy temp buffer to output buffer
    if(nSubdivisions > 0 && nSubdivisions % 2 != 0) {
      size_t n = nTriangles * 3;
      for(size_t i = 0; i < n; i++)
        connectivityList[i] = connectivityListTemp[i];
    }
  }

  // print progress
  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename DT, typename IT>
int ttk::Icosphere::translateIcosphere(DT *vertexCoords,
                                       IT *connectivityList,
                                       const size_t &icosphereIndex,
                                       const size_t &nVerticesPerIcosphere,
                                       const size_t &nTrianglesPerIcosphere,
                                       const DT *centers) const {
  size_t vertexCoordOffset = icosphereIndex * nVerticesPerIcosphere * 3;
  size_t connectivityListOffset = icosphereIndex * nTrianglesPerIcosphere * 3;
  size_t temp = icosphereIndex * 3;
  const DT &centerX = centers[temp++];
  const DT &centerY = centers[temp++];
  const DT &centerZ = centers[temp];

  // vertex coords
  for(size_t i = 0, limit = nVerticesPerIcosphere * 3; i < limit;) {
    vertexCoords[vertexCoordOffset++] = vertexCoords[i++] + centerX;
    vertexCoords[vertexCoordOffset++] = vertexCoords[i++] + centerY;
    vertexCoords[vertexCoordOffset++] = vertexCoords[i++] + centerZ;
  }

  // connectivity list
  size_t vertexIdOffset = icosphereIndex * nVerticesPerIcosphere;
  for(size_t i = 0, limit = nTrianglesPerIcosphere * 3; i < limit;) {
    connectivityList[connectivityListOffset++]
      = connectivityList[i++] + vertexIdOffset;
    connectivityList[connectivityListOffset++]
      = connectivityList[i++] + vertexIdOffset;
    connectivityList[connectivityListOffset++]
      = connectivityList[i++] + vertexIdOffset;
  }

  return 1;
}

template <typename DT, typename IT>
int ttk::Icosphere::computeIcospheres(
  // Output
  DT *vertexCoords,
  IT *connectivityList,

  // Input
  const size_t &nSpheres,
  const size_t &nSubdivisions,
  const DT &radius,
  const DT *centers,

  // Optional Output
  DT *normals) const {

  if(nSpheres < 1) {
    this->printWrn("Number of input points smaller than 1.");
    return 1;
  }

  // compute number of vertices and triangles for one ico sphere
  size_t nVerticesPerIcosphere, nTrianglesPerIcosphere;
  if(!this->computeNumberOfVerticesAndTriangles(
       nVerticesPerIcosphere, nTrianglesPerIcosphere, nSubdivisions))
    return 0;

  // compute ico sphere around origin
  if(!this->computeIcosphere(vertexCoords, connectivityList, nSubdivisions))
    return 0;

  // store normals if requested
  if(normals != nullptr) {
    ttk::Timer t;
    this->printMsg("Computing Normals", 0, 0, this->threadNumber_,
                   ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nSpheres; i++) {
      size_t n = nVerticesPerIcosphere * 3;
      size_t offset = i * n;
      for(size_t j = 0; j < n; j++)
        normals[offset++] = vertexCoords[j];
    }
    this->printMsg(
      "Computing Normals", 1, t.getElapsedTime(), this->threadNumber_);
  }

  // Translating Icospheres
  ttk::Timer timer;
  const std::string transMsg = "Moving " + std::to_string(nSpheres)
                               + " Icospheres (R: " + std::to_string(radius)
                               + ")";
  this->printMsg(
    transMsg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // apply radius
  if(radius != 1.0)
    for(size_t i = 0, j = nVerticesPerIcosphere * 3; i < j; i++) {
      vertexCoords[i] *= radius;
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 1; i < nSpheres; i++) {
    this->translateIcosphere(vertexCoords, connectivityList, i,
                             nVerticesPerIcosphere, nTrianglesPerIcosphere,
                             centers);
  }

  // translate first icosphere
  this->translateIcosphere(vertexCoords, connectivityList, 0,
                           nVerticesPerIcosphere, nTrianglesPerIcosphere,
                           centers);
  // print status
  this->printMsg(transMsg, 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
