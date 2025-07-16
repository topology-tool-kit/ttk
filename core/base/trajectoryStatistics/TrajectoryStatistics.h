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
#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#endif


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
                std::vector<std::vector<int>>    &trajVertexId,   //input
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
                std::vector<std::vector<double>> gradientNorms,
                std::vector<std::vector<double>> &merge,
                const triangulationType *triangulation);

    inline void setInputScalars(std::vector<void *> &is) {
      inputData_ = is;
    }
    
    inline void setFiltreX(double filtre) {
      filtreX_ = filtre;
    }

    inline void setFiltreY(double filtre) {
      filtreY_ = filtre;
    }

    inline void setCosCol(double filtre) {
        cosCol_ = filtre;
    }

    inline void setMaxRadus(double filtre){
        maxRadus_ = filtre;
    }

    inline void setMaxFrameDist(int filtre){
        maxFrameDist_ = filtre;
    }
// 1) Structure pour stocker une fusion i->j
    struct FuseRecord {
      int i, j;           // trajectoire i fusionnée vers trajectoire j
      int endFrame;       // frame de fin de i
      int startFrame;     // frame de début de j
      int finalContrib;
    };
    int correctTrajectory(
        std::vector<std::vector<int>>    &trajTime,
        std::vector<std::vector<double>> &coordsX,
        std::vector<std::vector<double>> &coordsY,
        std::vector<std::vector<double>> &merge,
        std::vector<std::vector<double>> &newTraj,
        std::vector<FuseRecord> &fuseRecords
    );


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

    #ifdef TTK_ENABLE_EIGEN
    int linearRegression(
        const std::vector<int> &T,   // times t_i
        const std::vector<double> &X,   // positions x_i
        const std::vector<double> &Y,   // positions y_i
        std::vector<double> &newTraj
    );
    #endif



    int computeMeanUnitDirectionLinear(
      const std::vector<std::vector<double>> &newTraj,
      std::vector<double> &meanDx,
      std::vector<double> &meanDy,
      std::vector<double> &meanDz
    ); 




    std::vector<void *> inputData_{};

    double filtreX_;
    double filtreY_;
    double cosCol_;
    double maxRadus_;
    int maxFrameDist_;

  }; // TrajectoryStatistics class

} // namespace ttk



#ifdef TTK_ENABLE_EIGEN
int ttk::TrajectoryStatistics::linearRegression(
  const std::vector<int> &T,   // times t_i
  const std::vector<double> &X,   // positions x_i
  const std::vector<double> &Y,   // positions y_i
  std::vector<double> &newTraj
) {
  const int n = (int)T.size();
  Eigen::MatrixXd M(n, 2);
  Eigen::VectorXd vx(n), vy(n);
  for(int i = 0; i < n; ++i) {
    M(i,0) = T[i];
    M(i,1) = 1.0;
    vx(i) = X[i];
    vy(i) = Y[i];
  }
  // (Mᵀ M) β = Mᵀ v  ⇒ β = [a; b]
  Eigen::Vector2d bxv = (M.transpose()*M).ldlt().solve(M.transpose()*vx);
  Eigen::Vector2d byv = (M.transpose()*M).ldlt().solve(M.transpose()*vy);
  newTraj.push_back(bxv[0]); // ax 
  newTraj.push_back(byv[0]); // ay ; 
  newTraj.push_back(bxv[1]); // bx
  newTraj.push_back(byv[1]); // by
  newTraj.push_back(-1); //futurFinalId

  return 1;
}
#endif




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


int ttk::TrajectoryStatistics::computeMeanUnitDirectionLinear(
  const std::vector<std::vector<double>> &newTraj,
  std::vector<double> &meanDx,
  std::vector<double> &meanDy,
  std::vector<double> &meanDz
)  {
  const size_t nTraj = newTraj.size();
  meanDx.assign(nTraj, 0.0);
  meanDy.assign(nTraj, 0.0);
  meanDz.assign(nTraj, 0.0);

  for(size_t i = 0; i < nTraj; ++i) {
    const auto &coef = newTraj[i]; // [ax, ay, bx, by]
    if(coef.size() < 2) continue;
    double ax = coef[0];
    double ay = coef[1];
    // Parametric direction vector (dx,dy,dz) = (ax, ay, 1)
    double dx = ax;
    double dy = ay;
    double dz = 1.0;
    double mag = std::sqrt(dx*dx + dy*dy + dz*dz);
    if(mag > 0.0) {
      meanDx[i] = dx / mag;
      meanDy[i] = dy / mag;
      meanDz[i] = dz / mag;
    }
  }

  return 1;
}

int ttk::TrajectoryStatistics::correctTrajectory(
    std::vector<std::vector<int>>    &trajTime,
    std::vector<std::vector<double>> &coordsX,
    std::vector<std::vector<double>> &coordsY,
    std::vector<std::vector<double>> &merge,
    std::vector<std::vector<double>> &newTraj,
    std::vector<FuseRecord> &fuseRecords
){
    const int numTraj = static_cast<int>(trajTime.size());
   
   #ifdef TTK_ENABLE_EIGEN   
    for (int i=0; i<numTraj; i++) {
        linearRegression(trajTime[i],coordsX[i],coordsY[i], newTraj[i]);
    }
    #endif
    
    this->printMsg("initially i have : " + std::to_string(numTraj) + " unique traj");

    std::vector<double> meanDx(numTraj);
    std::vector<double> meanDy(numTraj);
    std::vector<double> meanDz(numTraj);

    computeMeanUnitDirectionLinear(newTraj, meanDx, meanDy, meanDz);

    

    // === Étape A : collecte des merges directs ===
    
    fuseRecords.reserve(numTraj);

    std::vector<char> usedAsStart(numTraj, false), usedAsEnd(numTraj, false);
    
    this->printMsg("cosCol_ = " + std::to_string(cosCol_) + " maxRadus_ = " + std::to_string(maxRadus_));

    const double similarityThreshold = cosCol_;
    const double maxLinkDist2        = maxRadus_; // distance² maxi tolérée

    for(int i = 0; i < numTraj; ++i) {
      if(usedAsStart[i] || trajTime[i].empty()) continue;

      // géométrie et temps de fin de i
      const int endFrame = trajTime[i].back();

      double bestDot   = similarityThreshold;
      double bestDist2 = std::numeric_limits<double>::infinity();
      int    bestJ     = -1;

      // cherche le j qui maximise dot tout en respectant la distance
      for(int j = 0; j < numTraj; ++j) {
        if(usedAsEnd[j] || j == i || trajTime[j].empty()) continue;

        const int startFrame = trajTime[j].front();
        // contrainte temporelle
        if(startFrame <= endFrame || startFrame - endFrame > maxFrameDist_) continue;

        // similarité de direction
        double dot = meanDx[i]*meanDx[j]
                   + meanDy[i]*meanDy[j]
                   + meanDz[i]*meanDz[j];
        if(dot < bestDot) continue;

        // distance² au point projeté
        const auto &coefI = newTraj[i];
        const auto &coefJ = newTraj[j];
        const double xTh = coefI[0] * startFrame + coefI[2];
        const double yTh = coefI[1] * startFrame + coefI[3];
        const double zTh = static_cast<double>(startFrame);
        const double xJ  = coefJ[0] * startFrame + coefJ[2];
        const double yJ  = coefJ[1] * startFrame + coefJ[3];
        const double zJ  = static_cast<double>(startFrame);

        const double dx = xJ - xTh;
        const double dy = yJ - yTh;
        const double dz = zJ - zTh;
        const double dist2 = dx*dx + dy*dy + dz*dz;
        if(dist2 > maxLinkDist2) continue;

        // on garde si c'est mieux
        if (dist2 < bestDist2){
            bestDot   = dot;
            bestDist2 = dist2;
            bestJ     = static_cast<int>(j);
        }
      }

      // enregistrement si on a trouvé un match
      if(bestJ >= 0) {
        fuseRecords.push_back({ i,
                                bestJ,
                                trajTime[i].back(),
                                trajTime[bestJ].front(),
                                -1});
        usedAsStart[i] = true;
        usedAsEnd  [bestJ] = true;
      }
    }

    this->printMsg("j'ai reperé : " + std::to_string(fuseRecords.size()) + "fusion possible");

    merge.clear();
    merge.reserve(numTraj);
    std::vector<bool> used(fuseRecords.size(), false);
    int trajLost = 0;
    for (size_t idx1 = 0; idx1 < fuseRecords.size(); ++idx1) {
        if (used[idx1]) continue;
        auto &r1 = fuseRecords[idx1];
        int finalId = merge.size();
        std::vector<FuseRecord> finalTraj{r1};
        r1.finalContrib = finalId;
        newTraj[r1.i][4] = finalId;
        newTraj[r1.j][4] = finalId;
        used[idx1] = true;


        bool prepended = true;
        while (prepended) {
            prepended = false;
            for (size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
                if (used[idx2]) continue;
                auto &r2 = fuseRecords[idx2];
                if (r2.j == finalTraj.front().i) {
                    finalTraj.insert(finalTraj.begin(), r2); // insère au début
                    r2.finalContrib = finalId;
                    newTraj[r2.i][4] = finalId;
                    newTraj[r2.j][4] = finalId;
                    used[idx2] = true;
                    prepended = true;
                    break;
                }
            }
        }

        bool extended = true;
        while (extended) {
            extended = false;
            for (size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
              if (used[idx2]) continue;
              auto &r2 = fuseRecords[idx2];
              if (finalTraj.back().j == r2.i) {
                finalTraj.push_back(r2);
                r2.finalContrib = finalId;
                newTraj[r2.i][4] = finalId;
                newTraj[r2.j][4] = finalId;
                used[idx2] = true;
                extended = true;
                break;
              }
            }
        }
        
        int capacity = finalTraj.size()*2 + 2;  
        std::vector<int>    T;  T.reserve(capacity);
        std::vector<double> X;  X.reserve(capacity);
        std::vector<double> Y;  Y.reserve(capacity);

        for(const auto &r : finalTraj) {
          std::vector<int>    T2{trajTime[r.i].front(), r.endFrame};
          const auto &cI = newTraj[r.i];

          for (int t : T2){
            X.push_back(cI[0]*t + cI[2]);
            Y.push_back(cI[1]*t + cI[3]);
            T.push_back(t);
          }
        }
        
        FuseRecord &r = finalTraj.back();
        const auto &cJ = newTraj[r.j];
        X.push_back(cJ[0]*r.startFrame + cJ[2]);
        X.push_back(cJ[0]*trajTime[r.j].back() + cJ[2]);
        Y.push_back(cJ[1]*r.startFrame  + cJ[3]);
        Y.push_back(cJ[1]*trajTime[r.j].back() + cJ[3]);
        T.push_back(r.startFrame);
        T.push_back(trajTime[r.j].back());
        
        std::vector<double> lineCoef;
        linearRegression(T, X, Y, lineCoef);
        lineCoef[4] = trajTime[finalTraj[0].i].front();
        lineCoef.push_back(trajTime[r.j].back());
        
        double mag = std::sqrt(lineCoef[0]*lineCoef[0] + lineCoef[1]*lineCoef[1] + 1); 
        if ( (0.0 > lineCoef[0]/mag && lineCoef[0]/mag >= filtreX_) && (-filtreY_<lineCoef[1]/mag && lineCoef[1]/mag <= filtreY_)){
            merge.push_back(lineCoef);
        } else {
            trajLost++;
        }
    }
    this->printMsg("merge size avant traj non fus = " + std::to_string(merge.size()));
    this->printMsg("première condition a coupé : " +std::to_string(trajLost)); 
    for(int i = 0; i < numTraj; ++i) {
      if(!usedAsStart[i] && !usedAsEnd[i] && !trajTime[i].empty()) {
        //newTraj[i] == { ax, ay, bx, by } pour la trajectoire i
        if ( (0.0 > meanDx[i] && meanDx[i] >= filtreX_) && (-filtreY_<meanDy[i] && meanDy[i]<= filtreY_)){
            
            std::vector<double> lineCoef;
            lineCoef = newTraj[i];
            lineCoef[4] = trajTime[i].front();
            lineCoef.push_back(trajTime[i].back());
            newTraj[i][4] = merge.size();
            merge.push_back(lineCoef);

        } else {
            trajLost++;
        }
      }
    }
    this->printMsg("Au final j'ai : " + std::to_string(merge.size()) + " traj mais j'ai ai perdu : " + std::to_string(trajLost));
    this->printMsg("merge done");

    return 1;



}



template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::execute(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<int>>    &trajVertexId,
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
                std::vector<std::vector<double>> gradientNorms,
                std::vector<std::vector<double>> &finalTraj,
                const triangulationType *triangulation) {

    



    const int numTraj = static_cast<int>(trajTime.size());
    int nTraj = numTraj;
    const ttk::SimplexId numVertices = triangulation->getNumberOfVertices();
    const size_t numFrames = inputData_.size();
    const int  nPts = triangulation->getNumberOfVertices();

    
    const int numMerge = finalTraj.size();

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numMerge; ++i) {
        durations[i] = finalTraj[i][5] - finalTraj[i][4];
    }

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numMerge; ++i) {
        int dt = durations[i]; 
        VX[i] = (finalTraj[i][0] + finalTraj[i][2])/dt; // vx = (ax + bx) / (tf - t0)
        VY[i] = (finalTraj[i][0] + finalTraj[i][2])/dt;
    }

    // ####################### SURFACE ##########################

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
    



    this->printMsg("End base");
    return 1;
}



