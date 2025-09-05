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
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>


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
      triangulation->preconditionVertexStars();
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

    inline void setMinFrameDist(int filtre){
        minFrameDist_ = filtre;
    }

    inline void setSpatialScale(double filtre){
        spatialScale_ = filtre;
    }

    inline void setInterFrame(double filtre){
        interFrame_ = filtre;
    }
    
    inline void setConvertDur(bool filtre){
        convertDur_ = filtre;
    }

    inline void setMinVx(bool filtre){
        minVx_ = filtre;
    }

    inline void setCoordCratere(int filtre[2]){
        coordCratere_[0] = filtre[0];
        coordCratere_[1] = filtre[1];
    }

    inline void setThreshCratereAngle(double filtre){
        threshCratereAngle_ = filtre;
    }
    
    inline void setSurfaceMethod(int m){
        surfaceMethod_ = m;
    }

    inline void setMaxX(int filtre){
        maxX_ = filtre;
    }

    inline void setMaxY(int filtre){
        maxY_ = filtre;
    }

    inline void setMinY(int filtre){
        minY_ = filtre;
    }
    inline void setMinX(int filtre){
        minX_ = filtre;
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
    int bfsSegmentation(
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

    int computeSurfaceCellCount(const std::vector<ttk::SimplexId> &surfVertices,
                            const ttk::AbstractTriangulation *triangulation);

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

    #ifdef TTK_ENABLE_EIGEN
    template<class dataType>
    int randomWalkerSegment(
      const std::vector<ttk::SimplexId> &seed,        // ids des sommets "marqués"
      const std::vector<int> &seedLabel,              // label de chaque graine (0..K-1)
      const ttk::AbstractTriangulation *triangulation,      // maillage TTK (déjà initialisé)
      const dataType *intensities,                    // intensité par sommet
      const double beta,                               // paramètre des poids
      std::vector<int> &segmentation                   // [OUT] étiquette par sommet
    );
    #endif
    
    template<class dataType, class triangulationType>
    void collectNearMaxInSquare(
                           const triangulationType *tri,
                           const dataType *scalars,
                           const ttk::SimplexId centerId,
                           int square_size,
                           std::vector<ttk::SimplexId> &outIds
                           );

    template <class dataType, class triangulationType>
    int computeSurfacesBFS(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<int>>    &trajVertexId,
                std::vector<double>              &surfMin,
                std::vector<double>              &surfMax,
                std::vector<double>              &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                std::vector<ttk::SimplexId>      &excludedCriticalPoints,
                int                                frameSurf,
                double                             errSurf,
                std::vector<std::vector<double>>  gradientNorms,
                const triangulationType          *triangulation);


    template <class dataType, class triangulationType>
    int computeSurfacesRW(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<int>>    &trajVertexId,
                std::vector<double>              &surfMin,
                std::vector<double>              &surfMax,
                std::vector<double>              &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int                                frameSurf,
                double                             errSurf,
                std::vector<std::vector<double>>  gradientNorms,
                const triangulationType          *triangulation);

 

    std::vector<void *> inputData_{};

    double filtreX_;
    double filtreY_;
    double cosCol_;
    double maxRadus_;
    int maxFrameDist_;
    double spatialScale_;
    double interFrame_;
    bool convertDur_;
    double minVx_;
    int minFrameDist_;
    int coordCratere_[2];
    double threshCratereAngle_;
    int maxX_;
    int maxY_;
    int minY_;
    int minX_;
    int surfaceMethod_;
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

template<class dataType, class triangulationType>
void ttk::TrajectoryStatistics::collectNearMaxInSquare(
                           const triangulationType *tri,
                           const dataType *scalars,
                           const ttk::SimplexId centerId,
                           int square_size,
                           std::vector<ttk::SimplexId> &outIds
) {
  if(tri == nullptr || scalars == nullptr || square_size < 1)
    return;

  const auto nVerts = tri->getNumberOfVertices();
  if(nVerts <= 0 || centerId < 0 || centerId >= nVerts)
    return;


  // Coordonnées (arrondies) du centre
  float cxF=0.f, cyF=0.f, czF=0.f;
  tri->getVertexPoint(centerId, cxF, cyF, czF);
  const long long cx = llround(static_cast<double>(cxF));
  const long long cy = llround(static_cast<double>(cyF));

  const long long half = square_size / 2;

  // 1) Chercher le maximum dans le carré
  ttk::SimplexId maxId = -1;
  double maxVal = -std::numeric_limits<double>::infinity();

  for(ttk::SimplexId v = 0; v < nVerts; ++v) {
    float xF=0.f, yF=0.f, zF=0.f;
    tri->getVertexPoint(v, xF, yF, zF);
    const long long ix = llround(static_cast<double>(xF));
    const long long iy = llround(static_cast<double>(yF));

    if(std::llabs(ix - cx) > half || std::llabs(iy - cy) > half)
      continue;

    const double s = scalars[static_cast<size_t>(v)];
    if(s > maxVal) {
      maxVal = s;
      maxId = v;
    }
  }

  if(maxId < 0)
    return; // rien dans la fenêtre

  // 2) Pousser tous les sommets avec valeur dans [0.9*xmax, xmax]
  const double threshold = maxVal * 0.95;

  for(ttk::SimplexId v = 0; v < nVerts; ++v) {
    float xF=0.f, yF=0.f, zF=0.f;
    tri->getVertexPoint(v, xF, yF, zF);
    const long long ix = llround(static_cast<double>(xF));
    const long long iy = llround(static_cast<double>(yF));

    if(std::llabs(ix - cx) > half || std::llabs(iy - cy) > half)
      continue;

    const double s = scalars[static_cast<size_t>(v)];
    if(s >= threshold) {
      outIds.push_back(v); // inclut le maxId lui-même
    }
  }
}

int ttk::TrajectoryStatistics::computeSurfaceCellCount(
                            const std::vector<ttk::SimplexId> &surfVertices,
                            const ttk::AbstractTriangulation *triangulation) {
  std::unordered_set<ttk::SimplexId> cellIds;
  for(const ttk::SimplexId &v : surfVertices) {
    // Récupérer les cellules (étoiles) autour du sommet v
    const ttk::SimplexId starCount = triangulation->getVertexStarNumber(v);
    for(ttk::SimplexId k = 0; k < starCount; ++k) {
      ttk::SimplexId ttkCellId;
      triangulation->getVertexStar(v, k, ttkCellId);
      // Convertir en cellule VTK d'origine si applicable
      int vtkCellId;
      triangulation->getCellVTKID(ttkCellId, vtkCellId);
      cellIds.insert(vtkCellId);
    }
  }
  return static_cast<int>(cellIds.size());
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
    
    this->printMsg("N TRAJ TFF = " + std::to_string(numTraj));

    std::vector<double> meanDx(numTraj);
    std::vector<double> meanDy(numTraj);
    std::vector<double> meanDz(numTraj);

    computeMeanUnitDirectionLinear(newTraj, meanDx, meanDy, meanDz);

    
    fuseRecords.reserve(numTraj);

    std::vector<char> usedAsStart(numTraj, false), usedAsEnd(numTraj, false);

    const double similarityThreshold = cosCol_;
    const double maxLinkDist2        = maxRadus_; // distance² maxi tolérée

    for(int i = 0; i < numTraj; ++i) {
      if(usedAsStart[i] || trajTime[i].empty()) continue;

      const int endFrame = trajTime[i].back();

      double bestDot   = similarityThreshold;
      double bestDist2 = std::numeric_limits<double>::infinity();
      int    bestJ     = -1;

      // cherche le j qui maximise dot tout en respectant la distance
      for(int j = 0; j < numTraj; ++j) {
        if(usedAsEnd[j] || j == i || trajTime[j].empty()) continue;

        const int startFrame = trajTime[j].front();
        // contrainte temporelle
        if(startFrame - endFrame <= minFrameDist_  || startFrame - endFrame >= maxFrameDist_) continue;
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

    this->printMsg("N BOUT FUS = " + std::to_string(fuseRecords.size()));

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
        if (r1.j == 696) this->printMsg("being here");
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
       
        // ignoble a factoriser  

        double mag = std::sqrt(lineCoef[0]*lineCoef[0] + lineCoef[1]*lineCoef[1] + 1); 
        if ( (0.0 > lineCoef[0]/mag && lineCoef[0]/mag >= filtreX_) && (-filtreY_<lineCoef[1]/mag && lineCoef[1]/mag <= filtreY_) && std::abs(lineCoef[0]*spatialScale_*(1/interFrame_))>minVx_){
            if ((maxX_ != -1 && lineCoef[2] > maxX_) || (maxY_ != -1 && lineCoef[3] > maxY_) || (minY_ != -1 && lineCoef[3] < minY_ ) || (minX_ != -1 && lineCoef[2] < minX_))            
                for (int i=0; i<finalTraj.size(); i++){
                    newTraj[finalTraj[i].i][4] = -1;
                    newTraj[finalTraj[i].j][4] = -1;
                    trajLost++;
                }
            else {

                merge.push_back(lineCoef);
            }
        } else {
            trajLost++;
            for (int i=0; i<finalTraj.size(); i++){
                newTraj[finalTraj[i].i][4] = -1;
                newTraj[finalTraj[i].j][4] = -1;
            }
        }
    }
    this->printMsg("N TRAJ FUS " + std::to_string(merge.size()));
    this->printMsg("première condition a coupé : " +std::to_string(trajLost)); 
    
    for(int i = 0; i < numTraj; ++i) {
      if(!usedAsStart[i] && !usedAsEnd[i] && !trajTime[i].empty()) {
        //newTraj[i] == { ax, ay, bx, by } pour la trajectoire i
        if ( (0.0 > meanDx[i] && meanDx[i] >= filtreX_) && (-filtreY_<meanDy[i] && meanDy[i]<= filtreY_) && std::abs(newTraj[i][0]*spatialScale_*(1/interFrame_))>minVx_){
            std::vector<double> lineCoef;
            lineCoef = newTraj[i];
            lineCoef[4] = trajTime[i].front();
            lineCoef.push_back(trajTime[i].back());
            // produit scalaire entre traj et droite [cratère,pt de fin]

            double x_start = lineCoef[0]*lineCoef[4]+lineCoef[2];
            double y_start = lineCoef[1]*lineCoef[4]+lineCoef[3];
            double x_end = lineCoef[0]*lineCoef[5]+lineCoef[2];
            double y_end = lineCoef[1]*lineCoef[5]+lineCoef[3];
            double vx_traj = x_start - x_end;
            double vy_traj = y_start - y_end;
            double vx_crat = coordCratere_[0] - x_end;
            double vy_crat = coordCratere_[1] - y_end;
            double dot = vx_traj * vx_crat
                       + vy_traj * vy_crat;
            double traj_norm = std::sqrt(vx_traj*vx_traj + vy_traj*vy_traj);
            double crat_norm = std::sqrt(vx_crat*vx_crat + vy_crat*vy_crat);
            double angle = std::abs(dot/(traj_norm*crat_norm));
            if (maxX_ != -1 && lineCoef[2] > maxX_) {trajLost ++; continue;}
            if (maxY_ != -1 && lineCoef[3] > maxY_) {trajLost ++; continue;}
            if (minY_ != -1 && lineCoef[3] < minY_) {trajLost ++; continue;}
            if (minX_ != -1 && lineCoef[2] < minX_) {trajLost ++; continue;}

            if (angle>=threshCratereAngle_){
                newTraj[i][4] = merge.size();
                merge.push_back(lineCoef);
            }

        } else {
            trajLost++;
        }
      }
    }
    this->printMsg("N TRAJ TS = " + std::to_string(merge.size()));
    this->printMsg("N TRAJ SUPPR =" + std::to_string(trajLost));
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

    
    const int numMerge = finalTraj.size();

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    
    for(int i = 0; i < numMerge; ++i) {
        if (convertDur_)
            durations[i] = (finalTraj[i][5] - finalTraj[i][4])*interFrame_;
        else 
            durations[i] = finalTraj[i][5] - finalTraj[i][4];
    } 
   

    double conversion = spatialScale_*(1/interFrame_);
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numMerge; ++i) {
        int dt = durations[i]; 
        VX[i] = finalTraj[i][0]*conversion; // vx = ax 
        VY[i] = finalTraj[i][1]*conversion; // vy = ay
    }

    // ####################### SURFACE ##########################
    if(surfaceMethod_ == 0) {
      computeSurfacesBFS<dataType, triangulationType>(
        trajTime, trajVertexId,
        surfMin, surfMax, surfMoy,
        allVertexDebris, excludedCriticalPoints,
        frameSurf, errSurf, gradientNorms, triangulation);
    } else {
      computeSurfacesRW<dataType, triangulationType>(
        trajTime, trajVertexId,
        surfMin, surfMax, surfMoy,
        allVertexDebris, 
        frameSurf, errSurf, gradientNorms, triangulation);
    }
this->printMsg("End base");
    return 1;
}



template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::computeSurfacesBFS(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<int>>    &trajVertexId,
                std::vector<double>              &surfMin,
                std::vector<double>              &surfMax,
                std::vector<double>              &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                std::vector<ttk::SimplexId>      &excludedCriticalPoints,
                int                                frameSurf,
                double                             errSurf,
                std::vector<std::vector<double>>  gradientNorms,
                const triangulationType          *triangulation) {
// SURFACE (BFS)

    const size_t numFrames = inputData_.size();
    double maxVal = std::numeric_limits<double>::lowest();
    const int  nPts = triangulation->getNumberOfVertices();
    const int numTraj = static_cast<int>(trajTime.size());
    const ttk::SimplexId numVertices = triangulation->getNumberOfVertices();


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
    this->printMsg("entering");
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel num_threads(this->threadNumber_)
    {
      std::vector<char> visited(numVertices);
      std::vector<int> surfVertex;
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
                bfsSegmentation(vid, surfVertex, frameScalars, visited, local_min, errSurf, maxVal, gradientNorms[frame], triangulation);
                if(surfVertex.size() > 100) {
                    if(frame == frameSurf) {
                        excludedLocal[i] = vid;
                    }
                    surfVertex.clear();
                }

                if(frame == frameSurf) {
                    allVertexDebris[i] = surfVertex;
                }
                trajSurfaces[j] = computeSurfaceCellCount(surfVertex, triangulation);
            }

            auto [minIt, maxIt] = std::minmax_element(trajSurfaces.begin(), trajSurfaces.end());
            surfMin[i] = *minIt;
            surfMax[i] = *maxIt;
            long sum = 0;
            int count = 0;
            for(long v : trajSurfaces) {
              if(v != 0) {
                sum   += v;
                ++count;
              }
            }

            double mean = (count > 0) ? static_cast<double>(sum) / count : 0.0;

            if (mean == 0){
                surfMoy[i] = 0.25;
            } else
                surfMoy[i] = mean ;
            
      }
    } 
   
    this->printMsg("allVertexDebris size = " + std::to_string(allVertexDebris.size()));
    excludedCriticalPoints.clear();
    excludedCriticalPoints.reserve(numTraj);
    for(int i = 0; i < numTraj; ++i) {
        if(excludedLocal[i] != -1) {
            excludedCriticalPoints.push_back(excludedLocal[i]);
        }
    }
    


    return 0;
}


template<class dataType>
int ttk::TrajectoryStatistics::bfsSegmentation(
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

  // Calcul du seuil scalaire
  double coeff = (-1.0 * errSurf) / maxVal;
  const double scalarThreshold = local_min + (coeff * local_min + errSurf);

  // Coefficient pour tolérance sur sigma du gradient
  const double kSigma = 1.5;
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
      // Pas de voisin déjà accepté : on ne propage pas -> normalement
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
  }

  return anyAdded ? 1 : 0;
}


template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::computeSurfacesRW(
                std::vector<std::vector<int>>    &trajTime,
                std::vector<std::vector<int>>    &trajVertexId,
                std::vector<double>              &surfMin,
                std::vector<double>              &surfMax,
                std::vector<double>              &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int                                frameSurf,
                double                             errSurf,
                std::vector<std::vector<double>>  gradientNorms,
                const triangulationType          *triangulation) {
   
   
   const auto *frameScalars = static_cast<dataType *>(inputData_[frameSurf]);
   const int numTraj = static_cast<int>(trajTime.size());
   
   
   // 1) Collecte des seeds background 
    for(int i = 0; i < numTraj; i++) {
      for(size_t j = 0; j < trajVertexId[i].size(); j++) {
        if(trajTime[i][j] == frameSurf) {
          const ttk::SimplexId vid = trajVertexId[i][j];
          collectNearMaxInSquare(triangulation, frameScalars, vid, 30, allVertexDebris[i]);
        }
      }
    }

    // 2) Construire les seeds + labels (0 = background, i+1 = foreground de la trajectoire i)
    #ifdef TTK_ENABLE_EIGEN
    {
      std::vector<ttk::SimplexId> seed;  seed.reserve(1024);
      std::vector<int>            seedLabel; seedLabel.reserve(1024);

      const auto nVerts = triangulation->getNumberOfVertices();
      std::vector<char> isSeed(nVerts, 0); // pour dédupliquer

      // 2.a) BACKGROUND seeds 
      for(int i = 0; i < numTraj; i++) {
        for(const auto v : allVertexDebris[i]) {
          if(v >= 0 && v < nVerts && !isSeed[v]) {
            seed.push_back(v);
            seedLabel.push_back(0);
            isSeed[v] = 1;
          }
        }
      }

      // 2.b) FOREGROUND seeds = un seed par trajectoire à frameSurf (label = i+1)
      for(int i = 0; i < numTraj; i++) {
        for(size_t j = 0; j < trajVertexId[i].size(); j++) {
          if(trajTime[i][j] == frameSurf) {
            const auto v = static_cast<ttk::SimplexId>(trajVertexId[i][j]);
            if(v >= 0 && v < nVerts && !isSeed[v]) {
              seed.push_back(v);
              seedLabel.push_back(i + 1); 
              isSeed[v] = 1;
            }
          }
        }
      }

      // 3) Random Walker
      std::vector<int> segmentation; // sortie: label par sommet
      this->printMsg("RandomWalker: seeds=" + std::to_string(seed.size())
                     + ", beta(errSurf)=" + std::to_string(errSurf));

      //errSurf = beta
      const int rwStatus = randomWalkerSegment(
          seed, seedLabel, triangulation, frameScalars, static_cast<double>(errSurf), segmentation);

      if(rwStatus != 0) {
        this->printMsg("randomWalkerSegment failed with code " + std::to_string(rwStatus));
      } else {
        for(int i = 0; i < numTraj; i++) {
          allVertexDebris[i].clear();
        }
        for(ttk::SimplexId v = 0; v < static_cast<ttk::SimplexId>(segmentation.size()); v++) {
          const int lab = segmentation[v];
          if(lab > 0) { 
            const int trajIdx = lab - 1;
            if(trajIdx >= 0 && trajIdx < numTraj) {
              allVertexDebris[trajIdx].push_back(v);
            }
          }
        }
      }
    }
    #endif

    
    return 0;
}


#ifdef TTK_ENABLE_EIGEN
template<class dataType>
int ttk::TrajectoryStatistics::randomWalkerSegment(
  const std::vector<ttk::SimplexId> &seed,        // ids des sommets "marqués"
  const std::vector<int> &seedLabel,              // label de chaque graine (0..K-1)
  const ttk::AbstractTriangulation *triangulation,      // maillage TTK (déjà initialisé)
  const dataType *intensities,                    // intensité par sommet
  const double beta,                              // paramètre des poids
  std::vector<int> &segmentation                  // [OUT] étiquette par sommet
) {

  this->printMsg("RandomWalker: début de la fonction");

  if(!triangulation) {
    this->printMsg("ERROR : randomWalkerSegment: null triangulation.");
    return -1;
  }
  if(seed.size() != seedLabel.size()) {
    this->printMsg("ERROR : randomWalkerSegment: seed and seedLabel size mismatch.");
    return -2;
  }

  const auto nVerts = triangulation->getNumberOfVertices();
  if(nVerts <= 0) {
    this->printMsg("ERROR : randomWalkerSegment: empty triangulation.");
    return -3;
  }

  // --- 0) Préparation ---
  this->printMsg("RandomWalker: préparation des structures de données");
  segmentation.assign(nVerts, -1);

  std::vector<int> vertexLabel(nVerts, -1);
  int K = 0;
  for(size_t k = 0; k < seed.size(); ++k) {
    const auto v = seed[k];
    if(v < 0 || v >= nVerts) {
      this->printMsg("ERROR : randomWalkerSegment: seed vertex out of range.");
      return -4;
    }
    vertexLabel[v] = seedLabel[k];
    K = std::max(K, seedLabel[k] + 1);
  }
  if(K <= 0) {
    this->printMsg("ERROR : randomWalkerSegment: no labels found in seeds.");
    return -5;
  }

  std::vector<ttk::SimplexId> uIndex(nVerts, -1);
  std::vector<ttk::SimplexId> uVerts;
  uVerts.reserve(nVerts);
  for(ttk::SimplexId v = 0; v < nVerts; ++v) {
    if(vertexLabel[v] < 0) {
      uIndex[v] = static_cast<ttk::SimplexId>(uVerts.size());
      uVerts.push_back(v);
    }
  }
  const ttk::SimplexId nU = static_cast<ttk::SimplexId>(uVerts.size());
  if(nU == 0) {
    for(ttk::SimplexId v = 0; v < nVerts; ++v)
      segmentation[v] = vertexLabel[v];
    this->printMsg("RandomWalker: aucun nœud inconnu (tout est graine)");
    return 0;
  }

  // --- 1) Assemblage ---
  this->printMsg("RandomWalker: assemblage du Laplacien restreint et des RHS");

  using T = double;
  using Triplet = Eigen::Triplet<T>;
  std::vector<Eigen::VectorXd> rhs(K, Eigen::VectorXd::Zero(nU));
  std::vector<Triplet> L_triplets;

#ifdef TTK_ENABLE_OPENMP
  int nThreads = 1;
  #include <omp.h>
  nThreads = omp_get_max_threads();
  std::vector<std::vector<Triplet>> L_triplets_tls(static_cast<size_t>(nThreads));

  #pragma omp parallel for schedule(static)
  for(ttk::SimplexId ui = 0; ui < nU; ++ui) {
    const auto vi = uVerts[ui];
    T diag = 0.0;
    const int tid =
    #ifdef _OPENMP
      omp_get_thread_num();
    #else
      0;
    #endif
    auto &localTriplets = L_triplets_tls[static_cast<size_t>(tid)];

    const auto nNeigh = triangulation->getVertexNeighborNumber(vi);
    for(int ln = 0; ln < nNeigh; ++ln) {
      ttk::SimplexId vj{};
      triangulation->getVertexNeighbor(vi, ln, vj);

      const T gi = static_cast<T>(intensities[vi]);
      const T gj = static_cast<T>(intensities[vj]);
      const T diff = gi - gj;
      const T wij = std::exp(-beta * diff * diff);

      diag += wij;
      const auto uj = uIndex[vj];
      if(uj >= 0) {
        localTriplets.emplace_back(ui, uj, -wij);
      } else {
        const int lab = vertexLabel[vj];
        if(lab >= 0) {
          rhs[lab](ui) += wij;
        }
      }
    }
    localTriplets.emplace_back(ui, ui, diag);
  }

  for(auto &vec : L_triplets_tls) {
    L_triplets.insert(L_triplets.end(),
                      std::make_move_iterator(vec.begin()),
                      std::make_move_iterator(vec.end()));
  }
#else
  for(ttk::SimplexId ui = 0; ui < nU; ++ui) {
    const auto vi = uVerts[ui];
    T diag = 0.0;
    const auto nNeigh = triangulation->getVertexNeighborNumber(vi);
    for(int ln = 0; ln < nNeigh; ++ln) {
      ttk::SimplexId vj{};
      triangulation->getVertexNeighbor(vi, ln, vj);

      const T gi = static_cast<T>(intensities[vi]);
      const T gj = static_cast<T>(intensities[vj]);
      const T diff = gi - gj;
      const T wij = std::exp(-beta * diff * diff);

      diag += wij;
      const auto uj = uIndex[vj];
      if(uj >= 0) {
        L_triplets.emplace_back(ui, uj, -wij);
      } else {
        const int lab = vertexLabel[vj];
        if(lab >= 0) {
          rhs[lab](ui) += wij;
        }
      }
    }
    L_triplets.emplace_back(ui, ui, diag);
  }
#endif

  Eigen::SparseMatrix<T> L_U(nU, nU);
  L_U.setFromTriplets(L_triplets.begin(), L_triplets.end());
  L_U.makeCompressed();

  // --- 2) Factorisation ---
  this->printMsg("RandomWalker: factorisation du Laplacien");

  Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> llt;
  llt.compute(L_U);
  const bool useDirect = (llt.info() == Eigen::Success);

  Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<T>>
    cg;
  if(!useDirect) {
    cg.setMaxIterations(std::max<ttk::SimplexId>(2000, 5 * nU));
    cg.setTolerance(1e-10);
    cg.compute(L_U);
    if(cg.info() != Eigen::Success) {
      this->printMsg("ERROR : randomWalkerSegment: solver setup failed.");
      return -6;
    }
  }

  // --- 3) Résolution ---
  this->printMsg("RandomWalker: résolution des systèmes linéaires");

  std::vector<Eigen::VectorXd> X(K, Eigen::VectorXd::Zero(nU));
  for(int s = 0; s < K; ++s) {
    if(useDirect) {
      X[s] = llt.solve(rhs[s]);
      if(llt.info() != Eigen::Success) {
        this->printMsg("ERROR : randomWalkerSegment: LLT solve failed for label "
                       + std::to_string(s));
        return -7;
      }
    } else {
      X[s] = cg.solve(rhs[s]);
      if(cg.info() != Eigen::Success) {
        this->printMsg("ERROR : randomWalkerSegment: CG solve failed for label "
                       + std::to_string(s));
        return -8;
      }
    }
  }

  // --- 4) Attribution ---
  this->printMsg("RandomWalker: attribution des labels");

#ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for schedule(static)
  for(long long k = 0; k < static_cast<long long>(seed.size()); ++k) {
    const auto v = seed[static_cast<size_t>(k)];
    segmentation[v] = vertexLabel[v];
  }

  #pragma omp parallel for schedule(static)
  for(ttk::SimplexId ui = 0; ui < nU; ++ui) {
    int bestLab = 0;
    T bestVal = X[0](ui);
    for(int s = 1; s < K; ++s) {
      const T val = X[s](ui);
      if(val > bestVal) {
        bestVal = val;
        bestLab = s;
      }
    }
    segmentation[uVerts[ui]] = bestLab;
  }
#else
  for(const auto v : seed) {
    segmentation[v] = vertexLabel[v];
  }
  for(ttk::SimplexId ui = 0; ui < nU; ++ui) {
    int bestLab = 0;
    T bestVal = X[0](ui);
    for(int s = 1; s < K; ++s) {
      const T val = X[s](ui);
      if(val > bestVal) {
        bestVal = val;
        bestLab = s;
      }
    }
    segmentation[uVerts[ui]] = bestLab;
  }
#endif

  this->printMsg("RandomWalker: terminé avec succès");
  return 0;
}

#endif

