#include <TrajectoryStatistics.h>
#include <Triangulation.h>
#include <numeric>

ttk::TrajectoryStatistics::TrajectoryStatistics() {
  this->setDebugMsgPrefix("TrajectoryStatistics");
}

int ttk::TrajectoryStatistics::findSurface(
  ttk::SimplexId                         startId,
  std::vector<ttk::SimplexId>            &surfVertex,
  const std::vector<std::vector<double>> &vertexScalars,
  std::vector<char>                      &visited,
  const double                           local_min,
  int                                    frame,
  double                                 errSurf,
  double                                 maxVal,
  const ttk::AbstractTriangulation       *triangulation
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
                ttk::AbstractTriangulation *triangulation) {

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

