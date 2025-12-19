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
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>
//  #include <FTMTreePP.h>
#include <ExTreeM.h>
// #include <OrderDisambiguation.h>
#include <PathCompression.h>
#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>
// #include <vector>
// #include <algorithm>
// #include <cmath>
// #include <limits>

#endif


namespace ttk {
  class TrajectoryStatistics : virtual public Debug {

  public:
    TrajectoryStatistics();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      triangulation->preconditionVertexStars();
      return triangulation->preconditionVertexNeighbors();
    }


    inline void setInputScalars(std::vector<void *> &is) { inputData_ = is; }
    inline void setInstantPersistence(const std::vector<std::vector<double>> &P) { instantPers_ = P; }
    inline void setFiltreY(double v) { filtreY_ = v; }
    inline void setCosCol(double v) { cosCol_ = v; }
    inline void setMaxRadius(double v) { maxRadus_ = v; }
    inline void setMaxFrameDist(int v) { maxFrameDist_ = v; minFrameDist_ = -v;}
    inline void setSpatialScale(double v) { spatialScale_ = v; }
    inline void setInterFrame(double v) { interFrame_ = v; }
    inline void setConvertDur(bool v) { convertDur_ = v; }
    inline void setMinVx(double v) { minVx_ = v; }
    inline void setMaxVx(double v) { maxVx_ = v; }
	inline void setEnableFilteringMinVx(int v) { enableFilteringMinVx_ = v;}
	inline void setEnableFilteringTimeOrigin(int v) { enableFilteringTimeOrigin_ = v; }
	inline void SetEnableFilteringDuration(int v) { enableFilteringDuration_ = v; }
	inline void setEnableFilteringCosY(int v) { enableFilteringCosY_ = v; }
	inline void setDuraMin(int v) { duraMin_ = v; }
	inline void setXOrigin(int v){ xOrigin_ = v; }
	inline void setMinTimeOrigin(int v){ minTimeOrigin_ = v; }
	inline void setMaxX(int v){ maxX_ = v; }
    inline void setMaxY(int v){ maxY_ = v; }
    inline void setMinY(int v){ minY_ = v; }
    inline void setMinX(int v){ minX_ = v; }
    inline void setSurfaceMethod(int m){ surfaceMethod_ = m; }
	inline void setPersisThresh(double m){ persistenceThreshold_ = m; }
	inline void setMinSeg(std::vector<ttk::SimplexId> &m){ minSeg_ = &m; }
	inline void setSaddleSeg(std::vector<ttk::SimplexId> &m){ saddleSeg_ = &m; }

    struct FuseRecord {
      int i, j;           
      int endFrame;       
      int startFrame;     
      int finalContrib;
    };


    template <class dataType, class triangulationType>
    int execute(std::vector<std::vector<int>> &trajTime,         
                std::vector<std::vector<int>>    &trajVertexId,   
                std::vector<int> &durations,           
                std::vector<double> &VX,
                std::vector<double> &VY,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int frameSurf,
                double errSurf,
                std::vector<std::vector<double>> gradientNorms,
                std::vector<std::vector<double>> &merge,
  				std::vector<std::vector<double>> &distance,
                const triangulationType *triangulation);

    
    int correctTrajectory(
        std::vector<std::vector<int>>    &trajTime,
        std::vector<std::vector<double>> &coordsX,
        std::vector<std::vector<double>> &coordsY,
        std::vector<std::vector<double>> &merge,
        std::vector<std::vector<double>> &newTraj,
        std::vector<FuseRecord> &fuseRecords
    );

   protected: 



    #ifdef TTK_ENABLE_EIGEN
    int linearRegression(
        const std::vector<int> &T,   
        const std::vector<double> &X,  
        const std::vector<double> &Y, 
        std::vector<double> &newTraj
    );
    #endif



    int computeMeanUnitDirectionLinear(
      const std::vector<std::vector<double>> &newTraj,
      std::vector<double> &meanDx,
      std::vector<double> &meanDy,
      std::vector<double> &meanDz
    );



    template <class dataType, class triangulationType>
    int computeSurfacesBFS(
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


    template <class dataType, class triangulationType>
    int bfsSegmentation(
                ttk::SimplexId                         vertexId,
                std::vector<ttk::SimplexId>            &surfVertex,
                const dataType                         *frameScalars,
                std::vector<char>                      &visited,
                const double                           threshold,
                double                                 errSurf,
                double                                 maxVal,
                std::vector<double>                    &gradientNorm,
                const triangulationType                *triangulation
    );


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
                const triangulationType          *triangulation);

	template <class dataType, class triangulationType>
	int computeMergeTree(
               const ttk::SimplexId frameSurf,
               const triangulationType *triangulation,
               std::vector<std::vector<double>> &finalTraj,
               std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
               std::vector<double>              &surfMin,
               std::vector<double>              &surfMax,
               std::vector<double>              &surfMoy,
  			   std::vector<std::vector<double>> &distance
	); 


    template<class dataType, class triangulationType>
    void collectNearMaxInSquare(
                           const triangulationType *tri,
                           const dataType *scalars,
                           const ttk::SimplexId centerId,
                           int square_size,
                           std::vector<ttk::SimplexId> &outIds
    );

    #ifdef TTK_ENABLE_EIGEN
    template<class dataType>
    int randomWalkerSegment(
      const std::vector<ttk::SimplexId> &seed,        
      const std::vector<int> &seedLabel,              
      const ttk::AbstractTriangulation *triangulation,      
      const dataType *intensities,                    
      const double beta,                              
      std::vector<int> &segmentation                   
    );
    #endif


    template <class dataType, class triangulationType>
    int computeSurfacesPersistence(
      std::vector<std::vector<int>>    &trajTime,
      std::vector<std::vector<int>>    &trajVertexId,
      std::vector<double>              &surfMin,
      std::vector<double>              &surfMax,
      std::vector<double>              &surfMoy,
      std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
      int                                frameSurf,
      const triangulationType          *triangulation
    );
    

    int computeSurfaceCellCount(const std::vector<ttk::SimplexId> &surfVertices,
                            const ttk::AbstractTriangulation *triangulation
    );

	template <typename dataType, typename triangulationType>
	void cleanDarkSegmentInPlace(
								   std::vector<ttk::SimplexId> &segmentVerts,
                                   const dataType *scalars,
                                   const triangulationType *triangulation,
                                   const int otsuBins
	);


	template <class dataType>
	dataType otsuThresholdLocal(
							       const std::vector<ttk::SimplexId> &verts,
                                   const dataType *scalars,
                                   const int nbins
	);


    std::vector<void *> inputData_{};
    std::vector<std::vector<double>> instantPers_;

    double filtreY_;
    double cosCol_;
    double maxRadus_;
    int maxFrameDist_;
    double spatialScale_;
    double interFrame_;
    bool convertDur_;
    double minVx_;
	double maxVx_;
	int enableFilteringMinVx_;
	int enableFilteringTimeOrigin_;
	int enableFilteringCosY_;
	int enableFilteringDuration_;
	int duraMin_;
	int xOrigin_;
	int minTimeOrigin_;
	int minFrameDist_;
    int maxX_;
    int maxY_;
    int minY_;
    int minX_;
    int surfaceMethod_;
	double persistenceThreshold_;
	std::vector<ttk::SimplexId> *minSeg_;
	std::vector<ttk::SimplexId> *saddleSeg_;
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


/**
 * Collect vertices in a square window around a centerId whose
 * scalar value is within [0.95 * max_in_window, max_in_window].*
 */
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

  float cxF=0.f, cyF=0.f, czF=0.f;
  tri->getVertexPoint(centerId, cxF, cyF, czF);
  const long long cx = llround(static_cast<double>(cxF));
  const long long cy = llround(static_cast<double>(cyF));

  const long long half = square_size / 2;

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
    return; 

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
      outIds.push_back(v); 
    }
  }
}


/**
 * Count the number of unique VTK cell ids incident to the given
 * vertex set (via vertex-star traversal). 
 */
int ttk::TrajectoryStatistics::computeSurfaceCellCount(
                            const std::vector<ttk::SimplexId> &surfVertices,
                            const ttk::AbstractTriangulation *triangulation) {
  std::unordered_set<ttk::SimplexId> cellIds;
  for(const ttk::SimplexId &v : surfVertices) {
    const ttk::SimplexId starCount = triangulation->getVertexStarNumber(v);
    for(ttk::SimplexId k = 0; k < starCount; ++k) {
      ttk::SimplexId ttkCellId;
      triangulation->getVertexStar(v, k, ttkCellId);
      int vtkCellId;
      triangulation->getCellVTKID(ttkCellId, vtkCellId);
      cellIds.insert(vtkCellId);
    }
  }
  return static_cast<int>(cellIds.size());
}


/**
 * Compute per-trajectory unit direction vectors from linear
 * coefficients (ax, ay) 
 */
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

/**
 * Fuse and linearize trajectories into longer segments under
 * direction, temporal-gap and spatial-distance constraints.
 */
int ttk::TrajectoryStatistics::correctTrajectory(
    std::vector<std::vector<int>>    &trajTime,
    std::vector<std::vector<double>> &coordsX,
    std::vector<std::vector<double>> &coordsY,
    std::vector<std::vector<double>> &merge, // output Final Trajectory
    std::vector<std::vector<double>> &newTraj, // output linear Regression
    std::vector<FuseRecord> &fuseRecords
){
  const int numTraj = static_cast<int>(trajTime.size());

  auto projX = [&](const std::vector<double> &c, int t) -> double { return c[0] * t + c[2]; };
  auto projY = [&](const std::vector<double> &c, int t) -> double { return c[1] * t + c[3]; };

  auto dirDot = [&](int i, int j,
                    const std::vector<double> &meanDx,
                    const std::vector<double> &meanDy,
                    const std::vector<double> &meanDz) -> double {
    return meanDx[i] * meanDx[j] + meanDy[i] * meanDy[j] + meanDz[i] * meanDz[j];
  };

  auto temporalOk = [&](int startFrame, int endFrame) -> bool {
    return (startFrame - endFrame > minFrameDist_) && (startFrame - endFrame < maxFrameDist_);
  };

  auto dist2AtStartFrame = [&](const std::vector<double> &coefI,
                               const std::vector<double> &coefJ,
                               int startFrame) -> double {
    const double xTh = projX(coefI, startFrame);
    const double yTh = projY(coefI, startFrame);
    const double zTh = static_cast<double>(startFrame);
    const double xJ  = projX(coefJ, startFrame);
    const double yJ  = projY(coefJ, startFrame);
    const double zJ  = static_cast<double>(startFrame);
    const double dx = xJ - xTh, dy = yJ - yTh, dz = zJ - zTh;
    return dx * dx + dy * dy + dz * dz;
  };

  auto resetContribChain = [&](const std::vector<FuseRecord> &chain) {
    for(const auto &fr : chain) {
      newTraj[fr.i][4] = -1;
      newTraj[fr.j][4] = -1;
    }
  };

  auto violatesBBox = [&](const std::vector<double> &c) -> bool {
    // c = [ax, ay, bx, by, start, end]
    if(maxX_ != -1 && c[2] > maxX_) return true;
    if(maxY_ != -1 && c[3] > maxY_) return true;
    if(minY_ != -1 && c[3] < minY_) return true;
    if(minX_ != -1 && c[2] < minX_) return true;
    return false;
  };

  // Filtre sur l'inclinaison en Y : on travaille sur ny = composante Y d'un vecteur unité
  auto passInclinationYNy = [&](double ny) -> bool {
    if(enableFilteringCosY_ == 0) { return true; }
    return (-filtreY_ <= ny && ny <= filtreY_);
  };

  // Filtre sur la vitesse en X, à partir de la pente ax
  auto passSpeedXAx = [&](double ax) -> bool {
    const double vx_abs = ax * spatialScale_ * (1.0 / interFrame_);
    return (enableFilteringMinVx_ == 0.0) ? true : (vx_abs >= minVx_ && vx_abs <= maxVx_);
  };

  auto passDirSpeedNyAx = [&](double ny, double ax) -> bool {
    if(enableFilteringMinVx_ == 0.0 && filtreY_ == 1.0) {
      return true;
    }
    if(!passInclinationYNy(ny)) {
      return false;
    }
    return passSpeedXAx(ax);
  };

  // c = [ax, ay, bx, by, start, end]
  auto passDirSpeed = [&](const std::vector<double> &c) -> bool {
    const double mag = std::sqrt(c[0] * c[0] + c[1] * c[1] + 1.0);
    const double ny  = c[1] / mag; 
    const double ax  = c[0];
    return passDirSpeedNyAx(ny, ax);
  };

  auto passDura = [&](const std::vector<double> &c) -> bool {
  	if (enableFilteringDuration_ == 0) { return true; }
	const int start = c[4];
	const int end = c[5];
	return (duraMin_ <= std::abs(end - start));
  };

  auto passTimeOrigin = [&](const std::vector<double> &c) -> bool{
	if (enableFilteringTimeOrigin_ == 0) { return true; }
    const double ax     = c[0];
    const double bx     = c[2];
    const double eps = 1e-8;

    // Trajectoire quasi horizontale en Y
    if(std::abs(ax) < eps) {
      return true;
    }
    
	const double tCross = (xOrigin_ - bx) / ax;

    return (tCross >= minTimeOrigin_);  	
  };


  auto buildSamplesForChain = [&](const std::vector<FuseRecord> &finalTraj,
                                  std::vector<int> &T, std::vector<double> &X, std::vector<double> &Y) {
    int capacity = static_cast<int>(finalTraj.size()) * 2 + 2;
    T.reserve(capacity); X.reserve(capacity); Y.reserve(capacity);

    for(const auto &r : finalTraj) {
      std::vector<int> T2{trajTime[r.i].front(), r.endFrame};
      const auto &cI = newTraj[r.i];
      for(const int t : T2) {
        X.push_back(projX(cI, t));
        Y.push_back(projY(cI, t));
        T.push_back(t);
      }
    }
    const FuseRecord &r = finalTraj.back();
    const auto &cJ = newTraj[r.j];
    X.push_back(projX(cJ, r.startFrame));
    X.push_back(projX(cJ, trajTime[r.j].back()));
    Y.push_back(projY(cJ, r.startFrame));
    Y.push_back(projY(cJ, trajTime[r.j].back()));
    T.push_back(r.startFrame);
    T.push_back(trajTime[r.j].back());
  };

  auto fitLineCoefForChain = [&](const std::vector<FuseRecord> &finalTraj) -> std::vector<double> {
    std::vector<int>    T;
    std::vector<double> X, Y;
    buildSamplesForChain(finalTraj, T, X, Y);

    std::vector<double> lineCoef;
    linearRegression(T, X, Y, lineCoef); // lineCoef = [ax, ay, bx, by, futurFinalId=-1]
    lineCoef[4] = trajTime[finalTraj[0].i].front();               // start
    lineCoef.push_back(trajTime[finalTraj.back().j].back());      // end
    return lineCoef;
  };

#ifdef TTK_ENABLE_EIGEN
  for(int i = 0; i < numTraj; ++i) {
    linearRegression(trajTime[i], coordsX[i], coordsY[i], newTraj[i]);
  }
#endif

  std::vector<double> meanDx(numTraj), meanDy(numTraj), meanDz(numTraj);
  computeMeanUnitDirectionLinear(newTraj, meanDx, meanDy, meanDz);

  fuseRecords.reserve(numTraj);
  std::vector<char> usedAsStart(numTraj, false), usedAsEnd(numTraj, false);

  const double similarityThreshold = cosCol_;
  const double maxLinkDist2        = maxRadus_;

  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || trajTime[i].empty()) continue;

    const int endFrame = trajTime[i].back();

    double bestDot   = similarityThreshold;
    double bestDist2 = std::numeric_limits<double>::infinity();
    int    bestJ     = -1;

    for(int j = 0; j < numTraj; ++j) {
      if(usedAsEnd[j] || j == i || trajTime[j].empty()) continue;

      const int startFrame = trajTime[j].front();
      if(!temporalOk(startFrame, endFrame)) continue;

      const double dot = dirDot(i, j, meanDx, meanDy, meanDz);
      if(dot < bestDot) continue;

      const double dist2 = dist2AtStartFrame(newTraj[i], newTraj[j], startFrame);
      if(dist2 > maxLinkDist2) continue;

      if(dist2 < bestDist2) {
        bestDot   = dot;
        bestDist2 = dist2;
        bestJ     = j;
      }
    }

    if(bestJ >= 0) {
      fuseRecords.push_back({i, bestJ, trajTime[i].back(), trajTime[bestJ].front(), -1});
      usedAsStart[i]   = true;
      usedAsEnd  [bestJ] = true;
    }
  }

  merge.clear();
  merge.reserve(numTraj);

  std::vector<bool> used(fuseRecords.size(), false);

  for(size_t idx1 = 0; idx1 < fuseRecords.size(); ++idx1) {
    if(used[idx1]) continue;

    auto &r1 = fuseRecords[idx1];
    const int finalId = static_cast<int>(merge.size()); 

    std::vector<FuseRecord> finalTraj{r1};
    r1.finalContrib = finalId;
    newTraj[r1.i][4] = finalId;
    newTraj[r1.j][4] = finalId;
    used[idx1] = true;

    // prepend
    bool prepended = true;
    while(prepended) {
      prepended = false;
      for(size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
        if(used[idx2]) continue;
        auto &r2 = fuseRecords[idx2];
        if(r2.j == finalTraj.front().i) {
          finalTraj.insert(finalTraj.begin(), r2);
          r2.finalContrib = finalId;
          newTraj[r2.i][4] = finalId;
          newTraj[r2.j][4] = finalId;
          used[idx2] = true;
          prepended = true;
          break;
        }
      }
    }

    // extend
    bool extended = true;
    while(extended) {
      extended = false;
      for(size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
        if(used[idx2]) continue;
        auto &r2 = fuseRecords[idx2];
        if(finalTraj.back().j == r2.i) {
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

    std::vector<double> lineCoef = fitLineCoefForChain(finalTraj);

    if(passDirSpeed(lineCoef) && passTimeOrigin(lineCoef) && passDura(lineCoef)) {
      if(violatesBBox(lineCoef)) {
        resetContribChain(finalTraj);
      } else {
        merge.push_back(lineCoef);
      }
    } else {
      resetContribChain(finalTraj);
    }
  }


  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || usedAsEnd[i] || trajTime[i].empty()) continue;

    const double nyOrphan = meanDy[i];      
    const double axOrphan = newTraj[i][0]; 
    if(passDirSpeedNyAx(nyOrphan, axOrphan)) {

      std::vector<double> lineCoef = newTraj[i];
      lineCoef[4] = trajTime[i].front();
      lineCoef.push_back(trajTime[i].back());

      if(violatesBBox(lineCoef)) {
        continue;
      }
	  if(!passTimeOrigin(lineCoef) || !passDura(lineCoef)) { continue ; }

      newTraj[i][4] = static_cast<double>(merge.size());
      merge.push_back(lineCoef);
    }
  }


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
                int frameSurf,
                double errSurf,
                std::vector<std::vector<double>> gradientNorms,
                std::vector<std::vector<double>> &finalTraj,
  				std::vector<std::vector<double>> &distance,
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
   

    double conversion = spatialScale_*(1/interFrame_); // spatialScale_ en m
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numMerge; ++i) {
        VX[i] = finalTraj[i][0]*conversion; // vx = ax 
        VY[i] = finalTraj[i][1]*conversion; // vy = ay
    }

    if(surfaceMethod_ == 0) {
      computeSurfacesBFS<dataType, triangulationType>(
        trajTime, trajVertexId,
        surfMin, surfMax, surfMoy,
        allVertexDebris,
        frameSurf, errSurf, gradientNorms, triangulation);
    } if (surfaceMethod_ == 1) {
      computeSurfacesRW<dataType, triangulationType>(
        trajTime, trajVertexId,
        surfMin, surfMax, surfMoy,
        allVertexDebris, 
        frameSurf, errSurf, triangulation);
    } else if(surfaceMethod_ == 2) {
      computeSurfacesPersistence<dataType, triangulationType>(
        trajTime, trajVertexId,
        surfMin, surfMax, surfMoy,
        allVertexDebris,
        frameSurf, triangulation);
    } else if (surfaceMethod_ == 3) {
		computeMergeTree<dataType, triangulationType>(
				frameSurf,
				triangulation,
				finalTraj,
				allVertexDebris,
				surfMin, surfMax, surfMoy, distance);

	}

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
    int                                frameSurf,
    double                             errSurf,
    std::vector<std::vector<double>>  gradientNorms,
    const triangulationType          *triangulation) {

  // SURFACE (BFS)

  const size_t numFrames = inputData_.size();
  const int    nPts      = triangulation->getNumberOfVertices();
  const int    numTraj   = static_cast<int>(trajTime.size());
  const ttk::SimplexId numVertices = triangulation->getNumberOfVertices();

  if((int)surfMin.size() < numTraj) surfMin.resize(numTraj, 0.0);
  if((int)surfMax.size() < numTraj) surfMax.resize(numTraj, 0.0);
  if((int)surfMoy.size() < numTraj) surfMoy.resize(numTraj, 0.0);
  if((int)allVertexDebris.size() < numTraj) allVertexDebris.resize(numTraj);

  double maxVal = std::numeric_limits<double>::lowest();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(max : maxVal)
#endif
  for(size_t v = 0; v < numFrames; ++v) {
    auto *scalars = static_cast<dataType *>(inputData_[v]);
    const double localMax = *std::max_element(scalars, scalars + nPts);
    if(localMax > maxVal) maxVal = localMax;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<char>              visited(numVertices);
    std::vector<ttk::SimplexId>    surfVertex;
#pragma omp for schedule(dynamic)
    for(int i = 0; i < numTraj; ++i) {
#else
    std::vector<char>              visited(numVertices);
    std::vector<ttk::SimplexId>    surfVertex;
    for(int i = 0; i < numTraj; ++i) {
#endif

      const int trajSize = static_cast<int>(trajVertexId[i].size());
      std::vector<int> trajSurfaces(trajSize, 0);

      for(int j = 0; j < trajSize; ++j) {
        const int frame          = trajTime[i][j];
        const ttk::SimplexId vid = static_cast<ttk::SimplexId>(trajVertexId[i][j]);

        auto *frameScalars = static_cast<dataType *>(inputData_[frame]);
        const double local_min = frameScalars[vid];

        std::fill(visited.begin(), visited.end(), 0);
        surfVertex.clear();

		this->printMsg("frame = " + std::to_string(frame) + " frameSurf = " + std::to_string(frameSurf));
        bfsSegmentation(
          vid,                      // seed
          surfVertex,               // out vertices
          frameScalars,             // scalars at 'frame'
          visited,
          local_min,
          errSurf,
          maxVal,
          gradientNorms[frame],
          triangulation
        );

        if(surfVertex.size() > 100) { surfVertex.clear(); }
        if(frame == frameSurf) { allVertexDebris[i] = surfVertex; }

        trajSurfaces[j] = computeSurfaceCellCount(surfVertex, triangulation);
      }

      int    minVal = std::numeric_limits<int>::max();
      int    maxValS = 0;
      long   sum = 0;
      int    count = 0;

      for(const int s : trajSurfaces) {
        if(s > 0) {
          if(s < minVal) minVal = s;
          if(s > maxValS) maxValS = s;
          sum += s;
          ++count;
        }
      }

      if(count > 0) {
        surfMin[i] = static_cast<double>(minVal);
        surfMax[i] = static_cast<double>(maxValS);
        const double mean = static_cast<double>(sum) / static_cast<double>(count);
        surfMoy[i] = (mean == 0.0 ? 0.25 : mean); // conserve ta règle spéciale
      } else {
        surfMin[i] = 0.0;
        surfMax[i] = 0.0;
        surfMoy[i] = 0.25;
      }
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif

  return 0;
}
            


/**
 * Region grow (BFS) from a seed vertex using a scalar interval
 * and a local gradient consistency test
*/
template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::bfsSegmentation(
  ttk::SimplexId                         startId,
  std::vector<ttk::SimplexId>            &surfVertex,
  const dataType                         *frameScalars,
  std::vector<char>                      &visited,
  const double                           local_min,
  double                                 errSurf,
  double                                 maxVal,
  std::vector<double>                    &gradientNorm,
  const triangulationType                *triangulation
) {
  surfVertex.clear();

  double coeff = (-1.0 * errSurf) / maxVal;
  const double scalarThreshold = local_min + (coeff * local_min + errSurf);
  const double kSigma = 1.5;
  const double eps = 1e-6;

  std::vector<ttk::SimplexId> stack;
  stack.reserve(128);
  stack.push_back(startId);

  const size_t nPts = gradientNorm.size();
  std::vector<char> inSurf(nPts, 0);

  bool anyAdded = false;

  while(!stack.empty()) {
    const ttk::SimplexId vId = stack.back();
    stack.pop_back();

    if(visited[vId]) continue;
    visited[vId] = 1;

    // critical point always in 
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
      if(delta <= kSigma * sigmaGrad) gradAccepted = true;
    } else {
      // sigmaGrad ≈ 0 : homogeneous neighbor
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
    const triangulationType          *triangulation) {

  const auto *frameScalars = static_cast<dataType *>(inputData_[frameSurf]);
  const int numTraj = static_cast<int>(trajTime.size());

  if((int)surfMin.size() < numTraj) surfMin.resize(numTraj, 0.0);
  if((int)surfMax.size() < numTraj) surfMax.resize(numTraj, 0.0);
  if((int)surfMoy.size() < numTraj) surfMoy.resize(numTraj, 0.0);
  if((int)allVertexDebris.size() < numTraj) allVertexDebris.resize(numTraj);

  // 1)seeds background 
  for(int i = 0; i < numTraj; i++) {
    for(size_t j = 0; j < trajVertexId[i].size(); j++) {
      if(trajTime[i][j] == frameSurf) {
        const ttk::SimplexId vid = static_cast<ttk::SimplexId>(trajVertexId[i][j]);
        collectNearMaxInSquare(triangulation, frameScalars, vid, 30, allVertexDebris[i]);
      }
    }
  }

#ifdef TTK_ENABLE_EIGEN
  // 2) seeds + labels (0 = background, i+1 = foreground trajectory i)
  {
    std::vector<ttk::SimplexId> seed;     seed.reserve(1024);
    std::vector<int>            seedLabel; seedLabel.reserve(1024);

    const auto nVerts = triangulation->getNumberOfVertices();
    std::vector<char> isSeed(nVerts, 0); // déduplication

    // 2.a) BACKGROUND seeds (label 0)
    for(int i = 0; i < numTraj; i++) {
      for(const auto v : allVertexDebris[i]) {
        if(v >= 0 && v < nVerts && !isSeed[v]) {
          seed.push_back(v);
          seedLabel.push_back(0);
          isSeed[v] = 1;
        }
      }
    }

    // 2.b) FOREGROUND seeds
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
    std::vector<int> segmentation; // label for each vertex 
    this->printMsg("RandomWalker: seeds=" + std::to_string(seed.size())
                   + ", beta(errSurf)=" + std::to_string(errSurf));

    // errSurf = beta
    const int rwStatus = randomWalkerSegment(
      seed, seedLabel, triangulation, frameScalars,
      static_cast<double>(errSurf), segmentation);

    if(rwStatus != 0) {
      this->printMsg("randomWalkerSegment failed with code " + std::to_string(rwStatus));
    } else {
      for(int i = 0; i < numTraj; i++) {
        allVertexDebris[i].clear();
      }
      for(ttk::SimplexId v = 0; v < static_cast<ttk::SimplexId>(segmentation.size()); v++) {
        const int lab = segmentation[v];
        if(lab > 0) { // lab = i+1
          const int trajIdx = lab - 1;
          if(trajIdx >= 0 && trajIdx < numTraj) {
            allVertexDebris[trajIdx].push_back(v);
          }
        }
      }

      for(int i = 0; i < numTraj; ++i) {
        const int surfCells = computeSurfaceCellCount(allVertexDebris[i], triangulation);
        const double val = static_cast<double>(surfCells);
        surfMin[i] = val;
        surfMax[i] = val;
        surfMoy[i] = val;
      }
    }
  }
#else
  (void)errSurf; (void)trajTime; (void)trajVertexId; (void)frameSurf;
#endif

  return 0;
}

/**
 * Multi-label Random Walker segmentation from seeds
*/ 
#ifdef TTK_ENABLE_EIGEN
template<class dataType>
int ttk::TrajectoryStatistics::randomWalkerSegment(
  const std::vector<ttk::SimplexId> &seed,        // ids des sommets "marqués"
  const std::vector<int> &seedLabel,              // label de chaque graine (0..K-1)
  const ttk::AbstractTriangulation *triangulation, 
  const dataType *intensities,                    // intensité par sommet
  const double beta,                              // paramètre des poids
  std::vector<int> &segmentation                  // [OUT] label par sommet
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

template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::computeSurfacesPersistence(
  std::vector<std::vector<int>>                 &trajTime,
  std::vector<std::vector<int>>                 &trajVertexId,
  std::vector<double>                           &surfMin,
  std::vector<double>                           &surfMax,
  std::vector<double>                           &surfMoy,
  std::vector<std::vector<ttk::SimplexId>>      &allVertexDebris,
  int                                            frameSurf,
  const triangulationType                       *triangulation) {

  this->printMsg("Surface (Persistence/minima, multi-source) — frame "
                 + std::to_string(frameSurf));

  if(frameSurf < 0 || frameSurf >= static_cast<int>(inputData_.size())) {
    this->printMsg("computeSurfacesPersistence: invalid frameSurf");
    return -1;
  }
  if(!triangulation) {
    this->printMsg("computeSurfacesPersistence: null triangulation");
    return -2;
  }
  if(instantPers_.empty()) {
    this->printMsg("computeSurfacesPersistence: instantPers_ is empty");
    return -3;
  }

  const auto *frameScalars
    = static_cast<const dataType *>(inputData_[frameSurf]);
  const ttk::SimplexId nVerts
    = static_cast<ttk::SimplexId>(triangulation->getNumberOfVertices());

  struct Seed {
    int traj;               // index de trajectoire
    ttk::SimplexId v;       // sommet graine
    dataType fcrit;         // valeur au critique
    double pers;            // persistance instantanée (>= 0)
  };
  std::vector<Seed> seeds;
  seeds.reserve(trajTime.size());

  const int numTraj = static_cast<int>(trajTime.size());
  for(int i = 0; i < numTraj; ++i) {
    if(static_cast<size_t>(i) < allVertexDebris.size())
      allVertexDebris[i].clear();

    const auto &T = trajTime[i];
    const auto &V = trajVertexId[i];
    const auto &P = (i < static_cast<int>(instantPers_.size()))
                      ? instantPers_[i] : std::vector<double>{};

    for(size_t k = 0; k < T.size(); ++k) {
      if(T[k] == frameSurf) {
        const auto v = static_cast<ttk::SimplexId>(V[k]);
        if(v >= 0 && v < nVerts) {
          const dataType fcrit = frameScalars[v];
          double pers = (k < P.size() ? P[k] : 0.0);
          if(pers < 0.0) pers = 0.0;
          seeds.push_back({i, v, fcrit, pers});
        }
        break;
      }
    }
  }

  const size_t S = seeds.size();
  if(S == 0) {
    for(int i = 0; i < numTraj; ++i) {
      surfMin[i] = surfMax[i] = surfMoy[i] = 0.0;
    }
    this->printMsg("Surface (Persistence/minima) — no seeds at this frame");
    return 0;
  }

  //priorité / ties 
  std::vector<int> seedOrder(S);
  std::iota(seedOrder.begin(), seedOrder.end(), 0);
  std::stable_sort(seedOrder.begin(), seedOrder.end(),
                   [&](int a, int b){
                     if(seeds[a].fcrit != seeds[b].fcrit)
                       return seeds[a].fcrit < seeds[b].fcrit; // minima → plus bas d'abord
                     return seeds[a].traj < seeds[b].traj;
                   });
  std::vector<int> seedPrio(S, 0); // plus petit = plus prioritaire
  for(size_t rank = 0; rank < S; ++rank)
    seedPrio[seedOrder[rank]] = static_cast<int>(rank);

  //Selle associé à chaque seed
  std::vector<dataType> upper(S);
  for(size_t s = 0; s < S; ++s)
    upper[s] = static_cast<dataType>(seeds[s].fcrit + seeds[s].pers);

  struct QItem {
    dataType f;             // valeur du sommet candidat
    ttk::SimplexId v;       // sommet
    int s;                  // index du seed associé
    int prio;               // priorité du seed associé (plus petit gagne)
    bool operator<(QItem const &o) const {
      if(f != o.f) return f > o.f;          
      if(prio != o.prio) return prio > o.prio;
      return v > o.v;                       
    }
  };

  std::priority_queue<QItem> pq;
  std::vector<int> label(nVerts, -1);

  auto pushNeighbor = [&](ttk::SimplexId vj, int sIdx) {
    if(vj < 0 || vj >= nVerts) return;
    if(label[vj] != -1) return;
    const dataType fv = frameScalars[vj];
    if(fv <= upper[sIdx]) {
      pq.push(QItem{fv, vj, sIdx, seedPrio[sIdx]});
    }
  };

  for(size_t s = 0; s < S; ++s) {
    const auto v0 = seeds[s].v;
    if(v0 < 0 || v0 >= nVerts) continue;

    if(label[v0] == -1) 
      label[v0] = seeds[s].traj; // minimum surface = juste la seed

    const int nnei = triangulation->getVertexNeighborNumber(v0);
    for(int ln = 0; ln < nnei; ++ln) {
      ttk::SimplexId vj{-1};
      triangulation->getVertexNeighbor(v0, ln, vj);
      pushNeighbor(vj, static_cast<int>(s));
    }
  }

  while(!pq.empty()) {
    const auto it = pq.top(); pq.pop();
    const auto v  = it.v;
    const int s   = it.s;

    if(label[v] != -1) continue;          // déjà pris par un autre + prioritaire
    if(frameScalars[v] > upper[s]) continue;  // hors bande de la graine s (normalement pas possible)

    label[v] = seeds[s].traj;

    // propage aux voisins
    const int nnei = triangulation->getVertexNeighborNumber(v);
    for(int ln = 0; ln < nnei; ++ln) {
      ttk::SimplexId vj{-1};
      triangulation->getVertexNeighbor(v, ln, vj);
      pushNeighbor(vj, s);
    }
  }

  for(ttk::SimplexId v = 0; v < nVerts; ++v) {
    const int lab = label[v];
    if(lab >= 0 && lab < numTraj) {
      allVertexDebris[lab].push_back(v);
    }
  }

  for(int i = 0; i < numTraj; ++i) {
    const double surf = static_cast<double>(
      computeSurfaceCellCount(allVertexDebris[i], triangulation)
    );
    surfMin[i] = surfMax[i] = surfMoy[i] = surf;
  }

  this->printMsg("Surface (Persistence/minima, multi-source) — OK");
  return 0;
}


template <class dataType, class triangulationType>
int ttk::TrajectoryStatistics::computeMergeTree(
  const ttk::SimplexId frameSurf,
  const triangulationType *triangulation,
  std::vector<std::vector<double>> &finalTraj,
  std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
  std::vector<double>              &surfMin,
  std::vector<double>              &surfMax,
  std::vector<double>              &surfMoy,
  std::vector<std::vector<double>> &distance
) {

  const ttk::SimplexId nPixels = triangulation->getNumberOfVertices();
  const int nFrames = inputData_.size();
  std::vector<std::vector<double>> trajSurfaces(finalTraj.size());
  const auto nTraj = finalTraj.size();
  allVertexDebris.assign(nTraj, {});
  distance.assign(nTraj, {});
  std::vector<std::vector<double>> distanceFrame(nTraj);

  this->printMsg("Computing Merge Tree Segmentation");

  for(int frame = 0  ; frame < nFrames; frame++) {
	this->printMsg("Computing frame : " + std::to_string(frame));
    // -----------------------------------------------------------------------
    // 1) Persistence diagram : scalars + offsets (pattern TTK)
    // -----------------------------------------------------------------------
    auto *scalars = static_cast<dataType *>(inputData_[frame]);

    std::vector<ttk::SimplexId> pdOffsets(nPixels);
    ttk::preconditionOrderArray<dataType>(
      static_cast<size_t>(nPixels),
      scalars,
      pdOffsets.data(),
      this->threadNumber_);

    ttk::PersistenceDiagram persistenceDiagram;
    persistenceDiagram.setThreadNumber(this->threadNumber_);
    persistenceDiagram.setBackend(
      ttk::PersistenceDiagram::BACKEND::DISCRETE_MORSE_SANDWICH);
    persistenceDiagram.preconditionTriangulation(
      const_cast<triangulationType *>(triangulation));

    ttk::DiagramType diagram;

    const int statusPD = persistenceDiagram.execute(
      diagram,
      scalars,
      0, // H0
      pdOffsets.data(),
      const_cast<triangulationType *>(triangulation));

    if(statusPD != 0) {
      this->printErr("PersistenceDiagram::execute failed");
      return -1;
    }

    // -----------------------------------------------------------------------
    // 2) Extraction des sommets critiques (dimension 0, pers >= seuil)
    // -----------------------------------------------------------------------
    std::vector<ttk::SimplexId> criticalPoints;
    criticalPoints.reserve(diagram.size() * 2);

    for(const auto &pair : diagram) {
      if(pair.dim != 0)
        continue;

      const auto pers = pair.persistence();
      if(pers < this->persistenceThreshold_)
        continue;

      criticalPoints.push_back(pair.birth.id);
      criticalPoints.push_back(pair.death.id);
    }

    // -----------------------------------------------------------------------
    // 3) TopologicalSimplification
    // -----------------------------------------------------------------------
    const dataType *inputScalars = scalars;

    std::vector<dataType> outScalars(nPixels);
    std::copy(inputScalars, inputScalars + nPixels, outScalars.begin());

    std::vector<ttk::SimplexId> inputOffsets = pdOffsets;
    std::vector<ttk::SimplexId> offsets      = inputOffsets;

    ttk::TopologicalSimplification topoSimp;
    topoSimp.setThreadNumber(this->threadNumber_);
    topoSimp.setBackend(ttk::TopologicalSimplification::BACKEND::LTS);
    topoSimp.preconditionTriangulation(
      const_cast<triangulationType *>(triangulation));

    const bool addPerturbation = true;
    const ttk::SimplexId constraintNumber
      = static_cast<ttk::SimplexId>(criticalPoints.size());

    const ttk::DiagramType emptyDiagram;

    topoSimp.execute<dataType, triangulationType>(
      inputScalars,
      outScalars.data(),
      criticalPoints.empty() ? nullptr : criticalPoints.data(),
      inputOffsets.data(),
      offsets.data(),
      constraintNumber,
      addPerturbation,
      *const_cast<triangulationType *>(triangulation),
      emptyDiagram);

    // -----------------------------------------------------------------------
    // 4) Merge tree via ExTreeM
    // -----------------------------------------------------------------------

    std::vector<ttk::SimplexId> order(nPixels);
    ttk::preconditionOrderArray<dataType>(
      static_cast<size_t>(nPixels),
      outScalars.data(),
      order.data(),
      this->threadNumber_);

    std::vector<ttk::SimplexId> ascendingManifold(nPixels, -1);
    std::vector<ttk::SimplexId> descendingManifold(nPixels, -1);

    ttk::PathCompression pathComp;
    pathComp.setThreadNumber(this->threadNumber_);
    pathComp.setDebugLevel(this->debugLevel_);
    pathComp.setComputeSegmentation(true, true, false);

    ttk::PathCompression::OutputSegmentation om{
      ascendingManifold.data(),
      descendingManifold.data(),
      nullptr};

    {
      const int statusPC = pathComp.execute<triangulationType>(
        om,
        order.data(),
        *const_cast<triangulationType *>(triangulation));
      if(statusPC != 0) {
        this->printErr("PathCompression::execute failed");
        return -1;
      }
    }

    std::vector<ttk::SimplexId> segmentation(nPixels, -1);
    std::vector<char>           regionType(nPixels, 0);

    std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> persistencePairs;
    std::map<ttk::SimplexId, int> cpMap;
    std::vector<ttk::ExTreeM::Branch> branches;

    ttk::ExTreeM exTreeM;
    exTreeM.setThreadNumber(this->threadNumber_);
    exTreeM.setDebugLevel(this->debugLevel_);

    std::vector<ttk::SimplexId> orderJoin(order);

    // Inversion de l'ordre pour le Join Tree (parallélisée)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
    for(ttk::SimplexId i = 0; i < nPixels; ++i) {
#else
    for(ttk::SimplexId i = 0; i < nPixels; ++i) {
#endif
      orderJoin[i] = nPixels - orderJoin[i] - 1;
    }

    const auto treeType = ttk::ftm::TreeType::Join;

    int statusMT = exTreeM.computePairs<triangulationType>(
      persistencePairs,
      cpMap,
      branches,
      segmentation.data(),
      regionType.data(),
      ascendingManifold.data(),
      descendingManifold.data(),
      orderJoin.data(),
      triangulation,
      treeType);

    if(statusMT != 1) {
      this->printErr("ExTreeM::computePairs failed");
      return -1;
    }

    // -----------------------------------------------------------------------
    // 5) Association aux trajectoires
    // -----------------------------------------------------------------------


    std::vector<std::vector<ttk::SimplexId>> segmentId(
      *std::max_element(segmentation.begin(), segmentation.end()) + 1);

    for(size_t vId = 0; vId < segmentation.size(); vId++) {
      if(regionType[vId] == 0)
        segmentId[segmentation[vId]].push_back(static_cast<ttk::SimplexId>(vId));
    }

    std::vector<ttk::SimplexId> segMinVertex(segmentId.size(), -1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
    for(size_t s = 0; s < segmentId.size(); ++s) {
#else
    for(size_t s = 0; s < segmentId.size(); ++s) {
#endif
      const auto &verts = segmentId[s];
      if(verts.empty())
        continue;

      ttk::SimplexId bestV   = verts[0];
      auto           bestOrd = order[bestV];

      for(const auto v : verts) {
        const auto o = order[v];
        if(o > bestOrd) {
          bestOrd = o;
          bestV   = v;
        }
      }

      segMinVertex[s] = bestV;
    }

    for(size_t trajId = 0; trajId < nTraj; ++trajId) {
      const auto &traj = finalTraj[trajId];
      const double ax = traj[0];
      const double ay = traj[1];
      const double bx = traj[2];
      const double by = traj[3];
      const int    startF = static_cast<int>(traj[4]);
      const int    endF   = static_cast<int>(traj[5]);

      if(frame < startF || frame > endF)
        continue;

      const double x = ax * frame + bx;
      if(x < 0 || x > 384)
        continue;

      const double y = ay * frame + by;
      if(y < 0 || y > 224)
        continue;
      const ttk::SimplexId xi = std::lround(x);
      const ttk::SimplexId yi = std::lround(y);

      ttk::SimplexId vId = xi + yi * 384;

      if(vId < 0 || vId >= nPixels)
        continue;
	  
	  std::vector<char> segCleaned(segmentId.size(), 0);

	  if(regionType[vId] == 0) {
	    auto segId = segmentation[vId];

	    // --- POST-TRAITEMENT ---
		if(segId >= 0 && segId < (ttk::SimplexId)segmentId.size() && !segCleaned[segId] && segmentId[segId].size() > 8) 		 {
	  	  cleanDarkSegmentInPlace<dataType, triangulationType>(
	  	    segmentId[segId], scalars, triangulation, 64);
	  	  segCleaned[segId] = 1;
	    }

	    double surfVal = static_cast<double>(
	  	computeSurfaceCellCount(segmentId[segId], triangulation));
	    if(surfVal == 0) surfVal = 1;

	    trajSurfaces[trajId].push_back(surfVal);
	    distance[trajId].push_back(surfVal);
	    distanceFrame[trajId].push_back(frame);

	    if(frame == static_cast<int>(frameSurf)) {
	  	  allVertexDebris[trajId] = segmentId[segId];

	  	  (*saddleSeg_)[trajId]   = segMinVertex[segId];
	  	  (*minSeg_)[trajId]      = segId;
	    }
	  }

    } 
  }   



  // -------------------------------------------------------------------------
  // 6) Statistiques surfMin, surfMax, surfMoy par trajectoire (parallélisée)
  // -------------------------------------------------------------------------
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
  for(size_t trajId = 0; trajId < finalTraj.size(); ++trajId) {
#else
  for(size_t trajId = 0; trajId < finalTraj.size(); ++trajId) {
#endif
    const std::vector<double> &surfaces = trajSurfaces[trajId];

    double minVal = std::numeric_limits<double>::max();
    double maxVal = 0.0;
    double sum    = 0.0;
    int    count  = 0;

    for(const double s : surfaces) {
      if(s > 0) {
        if(s < minVal)
          minVal = s;
        if(s > maxVal)
          maxVal = s;
        sum += s;
        ++count;
      }
    }

    if(count > 0) {
      surfMin[trajId] = minVal;
      surfMax[trajId] = maxVal;
      const double mean = sum / static_cast<double>(count);
      surfMoy[trajId]   = mean;
    } else {
      surfMin[trajId] = 0.0;
      surfMax[trajId] = 0.0;
      surfMoy[trajId] = 0.0;
    }
  }

  return 0;
}


// ------------------------------------------------------------
// Helper: compute Otsu threshold on a set of scalar values
// Works for generic scalar types by discretizing into nbins.
// ------------------------------------------------------------
template <class dataType>
dataType ttk::TrajectoryStatistics::otsuThresholdLocal(
							       const std::vector<ttk::SimplexId> &verts,
                                   const dataType *scalars,
                                   const int nbins) {
  if(verts.empty())
    return dataType{0};

  dataType vmin = scalars[verts[0]];
  dataType vmax = scalars[verts[0]];
  for(const auto v : verts) {
    const auto s = scalars[v];
    if(s < vmin) vmin = s;
    if(s > vmax) vmax = s;
  }

  if(vmax <= vmin) // constant region
    return vmin;

  std::vector<double> hist(nbins, 0.0);

  const double minD = static_cast<double>(vmin);
  const double maxD = static_cast<double>(vmax);
  const double invRange = 1.0 / (maxD - minD);

  for(const auto v : verts) {
    const double s = static_cast<double>(scalars[v]);
    int b = static_cast<int>(std::floor((s - minD) * invRange * (nbins - 1)));
    b = std::max(0, std::min(nbins - 1, b));
    hist[b] += 1.0;
  }

  const double total = static_cast<double>(verts.size());
  for(auto &h : hist) h /= total;

  // cumulative sums
  std::vector<double> omega(nbins, 0.0); // weights
  std::vector<double> mu(nbins, 0.0);    // means
  omega[0] = hist[0];
  mu[0] = 0.0 * hist[0];
  for(int i = 1; i < nbins; ++i) {
    omega[i] = omega[i - 1] + hist[i];
    mu[i] = mu[i - 1] + static_cast<double>(i) * hist[i];
  }
  const double muT = mu[nbins - 1];

  int bestK = 0;
  double bestSigma = -1.0;

  for(int k = 0; k < nbins; ++k) {
    const double w0 = omega[k];
    const double w1 = 1.0 - w0;
    if(w0 <= 1e-12 || w1 <= 1e-12)
      continue;

    const double mu0 = mu[k] / w0;
    const double mu1 = (muT - mu[k]) / w1;

    const double sigmaB = w0 * w1 * (mu0 - mu1) * (mu0 - mu1);
    if(sigmaB > bestSigma) {
      bestSigma = sigmaB;
      bestK = k;
    }
  }

  // Map bin index back to scalar threshold
  const double t = minD + (static_cast<double>(bestK) / (nbins - 1)) * (maxD - minD);
  return static_cast<dataType>(t);
}

// ------------------------------------------------------------
// - compute local Otsu threshold T on this segment
// - keep only vertices with scalars[v] <= T (dark part)
// - keep only the largest connected component among kept vertices
// ------------------------------------------------------------
template <typename dataType, typename triangulationType>
void ttk::TrajectoryStatistics::cleanDarkSegmentInPlace(
								   std::vector<ttk::SimplexId> &segmentVerts,
                                   const dataType *scalars,
                                   const triangulationType *triangulation,
                                   const int otsuBins) {
  if(segmentVerts.size() < 2)
    return;

  // 1) Local adaptive threshold
  const dataType T = otsuThresholdLocal<dataType>(segmentVerts, scalars, otsuBins);

  // 2) Keep only dark vertices (<= T)
  std::vector<char> inSeg(triangulation->getNumberOfVertices(), 0);
  std::vector<ttk::SimplexId> kept;
  kept.reserve(segmentVerts.size());

  for(const auto v : segmentVerts) {
    if(scalars[v] <= T) {
      kept.push_back(v);
      inSeg[v] = 1;
    }
  }

  // If Otsu kept almost nothing, don't destroy the segment:
  // fall back to keep N darkest (robust fallback)
  const size_t segSize = segmentVerts.size();
  const size_t minKeep = std::max<size_t>(2, (size_t)std::ceil(0.10 * (double)segSize)); // 10%
  if(kept.size() < minKeep && segSize >= 2) {
  
    for(const auto v : kept) inSeg[v] = 0;
    kept.clear();

    std::vector<ttk::SimplexId> tmp = segmentVerts;
    std::sort(tmp.begin(), tmp.end(),
              [&](ttk::SimplexId a, ttk::SimplexId b) {
                return scalars[(ttk::SimplexId)a] < scalars[(ttk::SimplexId)b];
              });
    
    for(size_t i = 0; i < std::min(minKeep, tmp.size()); ++i) {
      kept.push_back(tmp[i]);
      inSeg[(ttk::SimplexId)tmp[i]] = 1;
    }
    
  }

  if(kept.empty())
    return;

  // 3) Connected components on kept vertices, keep the largest CC
  std::vector<char> visited(triangulation->getNumberOfVertices(), 0);

  std::vector<ttk::SimplexId> bestCC;
  bestCC.reserve(kept.size());

  std::queue<ttk::SimplexId> q;

  for(const auto seed : kept) {
    if(visited[seed])
      continue;

    // BFS
    std::vector<ttk::SimplexId> cc;
    cc.reserve(128);

    visited[seed] = 1;
    q.push(seed);

    while(!q.empty()) {
      const auto u = q.front();
      q.pop();
      cc.push_back(u);

      const auto deg = triangulation->getVertexNeighborNumber(u);
      for(ttk::SimplexId i = 0; i < deg; ++i) {
        ttk::SimplexId nb{};
        triangulation->getVertexNeighbor(u, i, nb);

        if(nb < 0) continue;
        if(!inSeg[nb]) continue;
        if(visited[nb]) continue;

        visited[nb] = 1;
        q.push(nb);
      }
    }

    if(cc.size() > bestCC.size())
      bestCC.swap(cc);
  }

  segmentVerts.swap(bestCC);
}

