/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TrajectoryStatistics
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %TrajectoryStatistics class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'TrajectoryStatistics'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The TrajectoryStatistics class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class TrajectoryStatistics : virtual public Debug {

  public:
    TrajectoryStatistics();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      return triangulation->preconditionVertexNeighbors();
    }
    
    int findSurface(
                    ttk::SimplexId                         vertexId,
                    std::vector<ttk::SimplexId>            &surfVertex,
                    const std::vector<std::vector<double>> &vertexScalars,
                    std::vector<char>                      &visited,
                    const double                           threshold,
                    int                                    frame,
                    double                                 errSurf,
                    double                                 maxVal,
                    const ttk::AbstractTriangulation    *triangulation
    ); 

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */

    template <class dataType, class triangulationType>
    int execute(std::vector<std::vector<int>> &trajTime,         //input
                std::vector<std::vector<double>> &trajX,         //input 
                std::vector<std::vector<double>> &trajY,         //input
                std::vector<std::vector<double>> &trajZ,         //input 
                std::vector<std::vector<int>>    &trajVertexId,   //input
                std::vector<std::vector<double>> &vertexScalars,
                std::vector<int> &startFrames,          //output
                std::vector<int> &endFrames,            //output
                std::vector<int> &durations,            //output
                std::vector<double> &VX,
                std::vector<double> &VY,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                std::vector<ttk::SimplexId> &excludedCriticalPoints,
                int frameSurf,
                double errSurf,
                const triangulationType *triangulation);

  }; // TrajectoryStatistics class

} // namespace ttk



int ttk::TrajectoryStatistics::findSurface(
  ttk::SimplexId                         startId,
  std::vector<ttk::SimplexId>            &surfVertex,
  const std::vector<std::vector<double>> &vertexScalars,
  std::vector<char>                      &visited,
  const double                           local_min,
  int                                    frame,
  double                                 errSurf,
  double                                 maxVal,
  const ttk::AbstractTriangulation      *triangulation
) {
  surfVertex.clear();
  std::vector<ttk::SimplexId> stack;
  stack.reserve(128);
  stack.push_back(startId);

  double coeff = (-1.0 * errSurf) / maxVal;
  const double threshold = local_min + (coeff * local_min + errSurf);
  bool anyAdded = false;

  // Parcours en profondeur (BFS) pour étendre la surface
  while(!stack.empty()) {
    auto vId = stack.back();
    stack.pop_back();

    if(visited[vId]) continue;
    visited[vId] = 1;

    double val = vertexScalars[vId][frame];
    if(val > threshold || val < local_min) continue; // local_min < val < threshold

    // Ce sommet est accepté dans la surface
    surfVertex.push_back(vId);
    anyAdded = true;

    const int nNbrs = triangulation->getVertexNeighborNumber(vId);
    for(int j = 0; j < nNbrs; ++j) {
      ttk::SimplexId nbr{-1};
      triangulation->getVertexNeighbor(vId, j, nbr);
      if(!visited[nbr]) {
        stack.push_back(nbr);
      }
    }
  }
  return anyAdded ? 1 : 0;
}

template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::execute(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<double>> &trajX,
                std::vector<std::vector<double>> &trajY,
                std::vector<std::vector<double>> &trajZ,
                std::vector<std::vector<int>>    &trajVertexId,
                std::vector<std::vector<double>> &vertexScalars,
                std::vector<int>                &startFrames,
                std::vector<int>                &endFrames,
                std::vector<int>                &durations,
                std::vector<double>             &VX,
                std::vector<double>             &VY,
                std::vector<double>             &surfMin,
                std::vector<double>             &surfMax,
                std::vector<double>             &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                std::vector<ttk::SimplexId>     &excludedCriticalPoints,
                int frameSurf,
                double errSurf,
                const triangulationType *triangulation) {

    const int numTraj = static_cast<int>(trajTime.size());
    const ttk::SimplexId numVertices = triangulation->getNumberOfVertices();

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numTraj; ++i) {
        startFrames[i] = trajTime[i][0];
        endFrames[i]   = trajTime[i].back();
        durations[i]   = endFrames[i] - startFrames[i];
    }

    const int z_translation = trajZ[0][1] - trajZ[0][0];
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numTraj; ++i) {
        double sommeVX = 0.0;
        double sommeVY = 0.0;
        const int pointCount = static_cast<int>(trajX[i].size());
        for(int j = 1; j < pointCount; ++j) {
            double dx = trajX[i][j]   - trajX[i][j-1];
            double dy = trajY[i][j]   - trajY[i][j-1];
            double dt = (trajZ[i][j]  - trajZ[i][j-1]) / static_cast<double>(z_translation);
            if(dt != 0.0) {
                sommeVX += dx / dt;
                sommeVY += dy / dt;
            }
        }
        double numPoints = static_cast<double>(pointCount) - 1.0;
        VX[i] = (numPoints != 0.0 ? sommeVX / numPoints : 0.0);
        VY[i] = (numPoints != 0.0 ? sommeVY / numPoints : 0.0);
    }

    double maxVal = std::numeric_limits<double>::lowest();
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_) reduction(max: maxVal)
    #endif
    for(size_t v = 0; v < vertexScalars.size(); ++v) {
        if(!vertexScalars[v].empty()) {
            double localMax = *std::max_element(vertexScalars[v].begin(), vertexScalars[v].end());
            if(localMax > maxVal) {
                maxVal = localMax;
            }
        }
    }
    this->printMsg("Max = " + std::to_string(maxVal));

    std::vector<ttk::SimplexId> excludedLocal(numTraj, -1);

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel num_threads(this->threadNumber_)
    {
      std::vector<char> visited(numVertices);
      std::vector<ttk::SimplexId> surfVertex;
      #pragma omp for schedule(dynamic)
      for(int i = 0; i < numTraj; ++i) {
    #else
      for(int i = 0; i < numTraj; ++i) {
    #endif

            const int trajSize = static_cast<int>(trajVertexId[i].size());
            std::vector<int> trajSurfaces(trajSize);
            for(int j = 0; j < trajSize; ++j) {
                const int frame          = trajTime[i][j];
                const ttk::SimplexId vid = static_cast<ttk::SimplexId>(trajVertexId[i][j]);
                const double local_min   = vertexScalars[vid][frame];

                std::fill(visited.begin(), visited.end(), 0);
                surfVertex.clear();

                findSurface(vid, surfVertex, vertexScalars, visited, local_min, frame, errSurf, maxVal, triangulation);

                if(surfVertex.size() > 500) {
                    surfVertex.clear();
                    if(frame == frameSurf) {
                        excludedLocal[i] = vid;
                    }
                }

                if(frame == frameSurf) {
                    allVertexDebris[i] = surfVertex;
                }
                trajSurfaces[j] = static_cast<int>(surfVertex.size());
            }

            auto [minIt, maxIt] = std::minmax_element(trajSurfaces.begin(), trajSurfaces.end());
            surfMin[i] = *minIt;
            surfMax[i] = *maxIt;
            long sum = std::accumulate(trajSurfaces.begin(), trajSurfaces.end(), 0l);
            surfMoy[i] = sum / static_cast<double>(trajSurfaces.size());
      }
    } 

    excludedCriticalPoints.clear();
    excludedCriticalPoints.reserve(numTraj);
    for(int i = 0; i < numTraj; ++i) {
        if(excludedLocal[i] != -1) {
            excludedCriticalPoints.push_back(excludedLocal[i]);
        }
    }

    this->printMsg("End base");
    return 1;
}
