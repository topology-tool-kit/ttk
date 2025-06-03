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
                const ttk::AbstractTriangulation    *triangulation
    );

    template <typename datatype>
    int computeGradientNorm(
      const std::vector<double> &x,
      const std::vector<double> &y,
      const std::vector<double> &z,
      std::vector<double> &gradientNorm,
      const datatype *scalarfield, // taille = npts
      const ttk::AbstractTriangulation *triangulation
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

    double val = frameScalars[vId];
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



template <typename datatype>
int ttk::TrajectoryStatistics::computeGradientNorm(
  const std::vector<double> &X,
  const std::vector<double> &Y,
  const std::vector<double> &Z,
  std::vector<double> &gradientNorm,
  const datatype *scalarField, // taille = npts
  const ttk::AbstractTriangulation *triangulation
) {
  const size_t nPts = X.size();

  #pragma omp parallel for num_threads(omp_get_max_threads())
  for(size_t v = 0; v < nPts; ++v) { 
    const int nNbrs = triangulation->getVertexNeighborNumber(v);
    if(nNbrs <= 0)
      continue;

    double gx = 0.0, gy = 0.0, gz = 0.0;
    double weightSum = 0.0;

    for(int j = 0; j < nNbrs; ++j) { // calcul gradient local (une direction)
      ttk::SimplexId nbr = -1;
      triangulation->getVertexNeighbor(v, j, nbr);

      double dx = X[nbr] - X[v];
      double dy = Y[nbr] - Y[v];
      double dz = Z[nbr] - Z[v]; // diff coord

      double distSq = dx * dx + dy * dy + dz * dz; //diff dot
      if(distSq == 0)
        continue;

      double w = 1.0 / distSq; // pondération
      double dv = static_cast<double>(scalarField[nbr]) - static_cast<double>(scalarField[v]); //diff scalar
    
      gx += w * dv * dx; // dérivée pour chaque direction
      gy += w * dv * dy; // avec somme sur tous les voisins
      gz += w * dv * dz;
      weightSum += w;
    }

    if(weightSum > 0.0) { // somme de tt les gradients locaux
      gx /= weightSum;
      gy /= weightSum; //normalisation
      gz /= weightSum;
      gradientNorm[v] = std::sqrt(gx * gx + gy * gy + gz * gz); //norme
      //if (v == 0)
        //this->printMsg("valeur en 0 calculée = " + std::to_string(gradientNorm[v]));
   } 
  }

  return 0;
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
        double neighborGrad = gradientNorm[nbr];
        double gradJump = std::abs(neighborGrad - gVal);

        // Seuillage sur le saut du gradient
        if(gradJump < gradientJumpThreshold) {
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

                    findSurface(vid, surfVertex, frameScalars, visited, local_min, errSurf, maxVal, triangulation);

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
    } else if (surfMethods == 2) {

      this->printMsg("GRADIANT METHODE");
      std::vector<ttk::SimplexId> excludedLocal(numTraj, -1);
      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel num_threads(this->threadNumber_)
      {
           std::vector<char> visited(numVertices);
           std::vector<ttk::SimplexId> surfVertex;
           std::vector<double> gradientNorm(nPts);
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
                std::fill(gradientNorm.begin(), gradientNorm.end(), 0);
                surfVertex.clear();
                //this->printMsg("frame = " + std::to_string(frame));
                computeGradientNorm(coordsX, coordsY, coordsZ, gradientNorm, frameScalars, triangulation);
                //this->printMsg("gradient norme 0 =" + std::to_string(gradientNorm[0]));
                int result;
                result = findSurfaceByGradient(vId, surfVertex, gradientNorm, visited,errSurf, triangulation);

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
