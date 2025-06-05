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



    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */

    template <class dataType, class triangulationType>
    int execute(std::vector<std::vector<int>> &trajTime,         //input
                std::vector<double> &coordsX,
                std::vector<double> &coordsY,
                std::vector<double> &coordsZ,
                std::vector<std::vector<int>>    &trajVertexId,   //input
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
                int surfMethods,
                std::vector<std::vector<double>> gradientNorms,
                const triangulationType *triangulation);
    inline void setInputScalars(std::vector<void *> &is) {
      inputData_ = is;
    }

   protected: 

    template <class dataType>
    int findSurface(
                ttk::SimplexId                         vertexId,
                std::vector<ttk::SimplexId>            &surfVertex,
                const dataType                         *frameScalars,
                std::vector<char>                      &visited,
                const double                           threshold,
                double                                 errSurf,
                double                                 maxVal,
                std::vector<double> &gradientNorm,
                const ttk::AbstractTriangulation    *triangulation
    );


    int findSurfaceByGradient(
      const ttk::SimplexId startId,
      std::vector<ttk::SimplexId> &surfVertex,
      std::vector<double> &gradientNorm,
      std::vector<char> &visited,
      double     errSurf,
      const ttk::AbstractTriangulation *triangulation
    ); 


    std::vector<void *> inputData_{};


  }; // TrajectoryStatistics class

} // namespace ttk


template<class dataType>
int ttk::TrajectoryStatistics::findSurface(
  ttk::SimplexId                         startId,
  std::vector<ttk::SimplexId>            &surfVertex,
  const dataType                         *frameScalars,
  std::vector<char>                      &visited,
  const double                           local_min,
  double                                 errSurf,
  double                                 maxVal,
  std::vector<double>                    &gradientNorm,
  const ttk::AbstractTriangulation       *triangulation
) {
  surfVertex.clear();

  // Calcul du seuil scalaire (cœur de la tâche)
  double coeff = (-1.0 * errSurf) / maxVal;
  const double scalarThreshold = local_min + (coeff * local_min + errSurf);

  // Coefficient pour tolérance sur sigma du gradient
  const double kSigma = 1.5;
  // Petite marge pour comparer des doubles
  const double eps = 1e-6;

  std::vector<ttk::SimplexId> stack;
  stack.reserve(128);
  stack.push_back(startId);

  // Marqueurs de visite et appartenance à la surface
  const size_t nPts = gradientNorm.size();
  std::vector<char> inSurf(nPts, 0);

  bool anyAdded = false;

  while(!stack.empty()) {
    const ttk::SimplexId vId = stack.back();
    stack.pop_back();

    if(visited[vId]) continue;
    visited[vId] = 1;

    // 1) Acceptation forcée pour le centre
    if(vId == startId) {
      surfVertex.push_back(vId);
      inSurf[vId] = 1;
      anyAdded = true;
      const int nNbrs0 = triangulation->getVertexNeighborNumber(vId);
      for(int j = 0; j < nNbrs0; ++j) {
        ttk::SimplexId nbr0{-1};
        triangulation->getVertexNeighbor(vId, j, nbr0);
        if(!visited[nbr0]) stack.push_back(nbr0);
      }
      continue;
    }

    // 2) Condition sur la valeur scalaire
    const double val = static_cast<double>(frameScalars[vId]);
    if(val >= local_min && val <= scalarThreshold) {
      surfVertex.push_back(vId);
      inSurf[vId] = 1;
      anyAdded = true;
      const int nNbrs = triangulation->getVertexNeighborNumber(vId);
      for(int j = 0; j < nNbrs; ++j) {
        ttk::SimplexId nbr{-1};
        triangulation->getVertexNeighbor(vId, j, nbr);
        if(!visited[nbr]) stack.push_back(nbr);
      }
      continue;
    }

    // 3) Condition sur le gradient, avec calcul sur voisins déjà acceptés
    const int nNbrs = triangulation->getVertexNeighborNumber(vId);
    double sumGrad = 0.0;
    std::vector<ttk::SimplexId> acceptedNbrs;
    acceptedNbrs.reserve(nNbrs);

    for(int j = 0; j < nNbrs; ++j) {
      ttk::SimplexId nbr{-1};
      triangulation->getVertexNeighbor(vId, j, nbr);
      if(inSurf[nbr]) {
        sumGrad += gradientNorm[nbr];
        acceptedNbrs.push_back(nbr);
      }
    }

    const size_t count = acceptedNbrs.size();
    if(count < 1) {
      this->printErr("IMPOSSIBLE");
      // Pas de voisin déjà accepté : on ne propage pas -> normalement
      // impossible
      continue;
    }

    const double meanGrad = sumGrad / static_cast<double>(count);
    double var = 0.0;
    for(const auto &nbr : acceptedNbrs) {
      double diff = gradientNorm[nbr] - meanGrad;
      var += diff * diff;
    }
    const double sigmaGrad = std::sqrt(var / static_cast<double>(count));
    const double currGrad = gradientNorm[vId];
    const double delta = std::abs(currGrad - meanGrad);

    bool gradAccepted = false;
    if(sigmaGrad > eps) {
      // tolérance normale
      if(delta <= kSigma * sigmaGrad) gradAccepted = true;
    } else {
      // sigmaGrad ≈ 0 : voisins très homogènes
      // on n'accepte le point que s'il est très proche
      if(delta <= eps) gradAccepted = true;
    }

    if(gradAccepted) {
      surfVertex.push_back(vId);
      inSurf[vId] = 1;
      anyAdded = true;
      for(int j = 0; j < nNbrs; ++j) {
        ttk::SimplexId nbr{-1};
        triangulation->getVertexNeighbor(vId, j, nbr);
        if(!visited[nbr]) stack.push_back(nbr);
      }
    }
    // Sinon on rejette ce sommet (pas de propagation)
  }

  return anyAdded ? 1 : 0;
}



int ttk::TrajectoryStatistics::findSurfaceByGradient(
  const ttk::SimplexId startId,
  std::vector<ttk::SimplexId> &surfVertex,
  std::vector<double> &gradientNorm,
  std::vector<char> &visited,
  double     errSurf,
  const ttk::AbstractTriangulation *triangulation
) {
  //this->printMsg("GRADIENT EDGE-BASED METHOD");
  constexpr size_t maxSurfaceSize = 100;
  double gradientJumpThreshold = errSurf; // à ajuster

  surfVertex.clear();
  std::vector<ttk::SimplexId> stack;
  stack.reserve(128);
  stack.push_back(startId);

  //const double refGrad = gradientNorm[startId];
  bool anyAdded = false;

  while(!stack.empty()) {
    auto vId = stack.back();
    stack.pop_back();

    if(visited[vId])
      continue;
    visited[vId] = 1;

    double gVal = gradientNorm[vId];

    // Ajout du point courant
    surfVertex.push_back(vId);
    anyAdded = true;

    if(surfVertex.size() > maxSurfaceSize) {
      surfVertex.clear(); // trop grand : rejet
      return 0;
    }

    const int nNbrs = triangulation->getVertexNeighborNumber(vId);
    for(int j = 0; j < nNbrs; ++j) {
      ttk::SimplexId nbr{-1};
      triangulation->getVertexNeighbor(vId, j, nbr);
      if(!visited[nbr]) {

        double ratio = gradientNorm[nbr] / (gVal + 1e-12);
        if(std::abs(1.0 - ratio) < gradientJumpThreshold) {
          stack.push_back(nbr);
        }

      }
    }
  }

  return anyAdded ? 1 : 0;
}



template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::execute(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<double> &coordsX,
                std::vector<double> &coordsY,
                std::vector<double> &coordsZ,
                std::vector<std::vector<int>>    &trajVertexId,
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
                int surfMethods,
                std::vector<std::vector<double>> gradientNorms,
                const triangulationType *triangulation) {

    const int numTraj = static_cast<int>(trajTime.size());
    const ttk::SimplexId numVertices = triangulation->getNumberOfVertices();
    const size_t numFrames = inputData_.size();
    const int  nPts = triangulation->getNumberOfVertices();


    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numTraj; ++i) {
        startFrames[i] = trajTime[i][0];
        endFrames[i]   = trajTime[i].back();
        durations[i]   = endFrames[i] - startFrames[i];
    }

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numTraj; ++i) {
        double sommeVX = 0.0;
        double sommeVY = 0.0;
        const int pointCount = static_cast<int>(trajVertexId[i].size());
        for(int j = 1; j < pointCount; ++j) {
            const int vId = trajVertexId[i][j];
            const int pvId = trajVertexId[i][j-1];
            double dx = coordsX[vId] - coordsX[pvId]; 
            double dy = coordsY[vId] - coordsY[pvId]; 



            sommeVX += dx;
            sommeVY += dy;
        }
        double numPoints = static_cast<double>(pointCount) - 1.0;
        VX[i] = (numPoints != 0.0 ? sommeVX / numPoints : 0.0);
        VY[i] = (numPoints != 0.0 ? sommeVY / numPoints : 0.0);
    }

    if (surfMethods == 1){
        double maxVal = std::numeric_limits<double>::lowest();
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(this->threadNumber_) reduction(max: maxVal)
        #endif
        for(size_t v = 0; v < numFrames; ++v) {
            auto *scalars = static_cast<dataType*>(inputData_[v]);
            double localMax = *std::max_element(scalars, scalars + nPts);
            if(localMax > maxVal) {
                maxVal = localMax;
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
                    auto *frameScalars = static_cast<dataType*>(inputData_[frame]);
                    const double local_min = frameScalars[vid];

                    std::fill(visited.begin(), visited.end(), 0);
                    surfVertex.clear();

                    findSurface(vid, surfVertex, frameScalars, visited, local_min, errSurf, maxVal, gradientNorms[frame], triangulation);

                    if(surfVertex.size() > 100) {
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
    } else if (surfMethods == 2) {

      this->printMsg("GRADIANT METHODE");
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
             for (int j=0; j<trajSize; j++){
                const int frame = trajTime[i][j];
                const ttk::SimplexId vId = static_cast<ttk::SimplexId>(trajVertexId[i][j]);
                auto *frameScalars = static_cast<dataType*>(inputData_[frame]);
                std::fill(visited.begin(), visited.end(), 0);
                surfVertex.clear();

                //this->printMsg("gradient norme 0 =" + std::to_string(gradientNorm[0]));
                int result;
                result = findSurfaceByGradient(vId, surfVertex, gradientNorms[frame], visited,errSurf, triangulation);

                if (result == 0){
                    if (frame == frameSurf){
                        excludedLocal[i] = vId;
                    }
                }
                
                if (frame == frameSurf) {
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
    }

    this->printMsg("End base");
    return 1;
}
