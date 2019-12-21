/// \ingroup base
/// \class ttk::HelloWorld
/// \author Your Name Here <Your Email Address Here>
/// \date The Data Here.
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  class HelloWorld : virtual public Debug {

  public:
    HelloWorld() {
      this->setDebugMsgPrefix(
        "HelloWorld"); // inherited from Debug: prefix will be printed at the
                       // beginning of every msg
    };
    ~HelloWorld(){};

    template <class idType>
    int computeBoundingBox(float *boundingBoxPointCoordinates,
                           idType *boundingBoxConnectivityList,
                           ttk::Triangulation *triangulation,
                           const float &scale) const;

  private:
    int preconditionTriangulation(ttk::Triangulation *triangulation) const;
  };
} // namespace ttk

int ttk::HelloWorld::preconditionTriangulation(
  ttk::Triangulation *triangulation) const {
  if(!triangulation->hasPreconditionedVertexNeighbors())
    return triangulation->preconditionVertexNeighbors();
  return 0;
};

template <class idType>
int ttk::HelloWorld::computeBoundingBox(float *boundingBoxPointCoordinates,
                                        idType *boundingBoxConnectivityList,
                                        ttk::Triangulation *triangulation,
                                        const float &scale) const {

  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator

// check general input parameter validity only if TTK_ENABLE_KAMIKAZE is
// disabled
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation || !boundingBoxPointCoordinates
     || !boundingBoxConnectivityList) {
    this->printErr(
      "Invalid triangulation, coordinates, or connectivity list pointer.");
    return 0; // return failure
  }
#endif

  // check algorithm specific input parameter validity
  if(triangulation->isEmpty() || triangulation->getNumberOfVertices() < 1) {
    this->printErr(
      "Unable to compute bounding box for triangulation without vertices.");
    return 0; // return failure
  }

  // ensure that triangulation is preconditioned for the following operations
  this->preconditionTriangulation(triangulation);

  // start a global timer AFTER preconditioning
  ttk::Timer globalTimer;

  // variables that will store BB extent
  float xMin, xMax, yMin, yMax, zMin, zMax = 0;

  // -------------------------------------------------------------------------
  // compute bounding box extent of triangulation
  // -------------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer timer;

    // print the progress of the current subprocedure (at the beginning 0)
    this->printMsg("Compute BB extent",
                   0 // progress form 0-1
    );

    // variables that will store the xyz coordinates of the current vertex
    float x, y, z = 0;

    // nVertices has to > 0 otherwise the triangulation would be empty
    size_t nVertices = triangulation->getNumberOfVertices();

    // initialize BB extent with first vertex
    triangulation->getVertexPoint(0, x, y, z);
    xMin = x;
    xMax = x;
    yMin = y;
    yMax = y;
    zMin = z;
    zMax = z;

    // compute actual extent
    for(size_t i = 1; i < nVertices; i++) {
      triangulation->getVertexPoint(i, x, y, z);

      xMin = std::min(xMin, x);
      yMin = std::min(yMin, y);
      zMin = std::min(zMin, z);

      xMax = std::max(xMax, x);
      yMax = std::max(yMax, y);
      zMax = std::max(zMax, z);
    }

    // print the progress of the current subprocedure with elapsed time
    this->printMsg(
      "Compute BB extent", 1, timer.getElapsedTime(), // progress, time
      ttk::debug::LineMode::REPLACE // replace last line of output stream
    );

    // in case of detailed reporting print extent to stream as a table
    this->printMsg({{"xBounds", "[" + std::to_string(xMin) + ", "
                                  + std::to_string(xMax) + "]"},
                    {"yBounds", "[" + std::to_string(yMin) + ", "
                                  + std::to_string(yMax) + "]"},
                    {"zBounds", "[" + std::to_string(zMin) + ", "
                                  + std::to_string(zMax) + "]"}},
                   ttk::debug::Priority::DETAIL);
  }

  // -------------------------------------------------------------------------
  // Compute corner point coordinates
  // -------------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer timer;

    // print the progress of the current subprocedure (at the beginning 0)
    this->printMsg("Compute BB corner points",
                   0 // progress form 0-1
    );

    // set corner points coordinates
    size_t index = 0;

    // Bottom 4 vertices
    // 0. vertex
    boundingBoxPointCoordinates[index++] = xMin;
    boundingBoxPointCoordinates[index++] = yMin;
    boundingBoxPointCoordinates[index++] = zMin;

    // 1. vertex
    boundingBoxPointCoordinates[index++] = xMax;
    boundingBoxPointCoordinates[index++] = yMin;
    boundingBoxPointCoordinates[index++] = zMin;

    // 2. vertex
    boundingBoxPointCoordinates[index++] = xMin;
    boundingBoxPointCoordinates[index++] = yMax;
    boundingBoxPointCoordinates[index++] = zMin;

    // 3. vertex
    boundingBoxPointCoordinates[index++] = xMax;
    boundingBoxPointCoordinates[index++] = yMax;
    boundingBoxPointCoordinates[index++] = zMin;

    // Top 4 vertices
    // 4. vertex
    boundingBoxPointCoordinates[index++] = xMin;
    boundingBoxPointCoordinates[index++] = yMin;
    boundingBoxPointCoordinates[index++] = zMax;

    // 5. vertex
    boundingBoxPointCoordinates[index++] = xMax;
    boundingBoxPointCoordinates[index++] = yMin;
    boundingBoxPointCoordinates[index++] = zMax;

    // 6. vertex
    boundingBoxPointCoordinates[index++] = xMin;
    boundingBoxPointCoordinates[index++] = yMax;
    boundingBoxPointCoordinates[index++] = zMax;

    // 7. vertex
    boundingBoxPointCoordinates[index++] = xMax;
    boundingBoxPointCoordinates[index++] = yMax;
    boundingBoxPointCoordinates[index++] = zMax;

// Multiply coordinates by scale (assumes object center is [0,0,0])
// Since each computation is independent we can use a parallel for loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < 24; i++)
      boundingBoxPointCoordinates[i] *= scale;

    // print the progress of the current subprocedure with elapsed time
    this->printMsg(
      "Compute BB corner points", 1, timer.getElapsedTime(), // progress, time
      ttk::debug::LineMode::REPLACE // replace last line of output stream
    );
  }

  // -------------------------------------------------------------------------
  // update connectivity list
  // -------------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer timer;

    // print the progress of the current subprocedure (at the beginning 0)
    this->printMsg("Compute BB connectivity list",
                   0 // progress form 0-1
    );

    size_t index = 0;

    // Voxel
    boundingBoxConnectivityList[index++] = 8; // #vertices that constitute cell
    boundingBoxConnectivityList[index++] = 0; // vertex index
    boundingBoxConnectivityList[index++] = 1; // vertex index
    boundingBoxConnectivityList[index++] = 2; // vertex index
    boundingBoxConnectivityList[index++] = 3; // vertex index
    boundingBoxConnectivityList[index++] = 4; // vertex index
    boundingBoxConnectivityList[index++] = 5; // vertex index
    boundingBoxConnectivityList[index++] = 6; // vertex index
    boundingBoxConnectivityList[index++] = 7; // vertex index

    // print the progress of the current subprocedure with elapsed time
    this->printMsg(
      "Compute BB connectivity list", 1,
      timer.getElapsedTime(), // progress, time
      ttk::debug::LineMode::REPLACE // replace last line of output stream
    );
  }

  // -------------------------------------------------------------------------
  // print global performance
  // -------------------------------------------------------------------------
  {
    this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
    this->printMsg(
      "Complete", 1, globalTimer.getElapsedTime() // global progress, time
    );
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1; // return success
}
