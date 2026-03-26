/// \ingroup base
/// \class ttk::DebrisTracer
/// \author Théophane Loloum <theophane.loloum@gmail.com> 
/// \date March 2026
///
/// This module defines the %DebrisTracer class that linearizes, chains, and
/// extracts surface statistics from tracked debris trajectories.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>
#include <ExTreeM.h>
#include <PathCompression.h>
#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#endif


namespace ttk {
  class DebrisTracer : virtual public Debug {

  public:
    DebrisTracer();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      return triangulation->preconditionVertexStars();
    }


    inline void setInputScalars(std::vector<void *> &inputScalars) { inputData_ = inputScalars; }
    inline void setInstantPersistence(const std::vector<std::vector<double>> &P) { instantPers_ = P; }
    inline void setFiltreY(double v) { filtreY_ = v; }
    inline void setCosCol(double v) { cosCol_ = v; }
    inline void setMaxRadius(double v) { maxRadius_ = v; }
    inline void setMaxFrameDist(int v) { maxFrameDist_ = v; minFrameDist_ = -v;}
    inline void setSpatialScale(double v) { spatialScale_ = v; }
    inline void setInterFrame(double v) { interFrame_ = v; }
    inline void setConvertDur(bool v) { convertDur_ = v; }
    inline void setMinVx(double v) { minVx_ = v; }
    inline void setMaxVx(double v) { maxVx_ = v; }
	inline void setEnableFilteringMinVx(bool v) { enableFilteringMinVx_ = v;}
	inline void setEnableFilteringTimeOrigin(bool v) { enableFilteringTimeOrigin_ = v; }
	inline void setEnableFilteringDuration(bool v) { enableFilteringDuration_ = v; }
	inline void setEnableFilteringCosY(bool v) { enableFilteringCosY_ = v; }
	inline void setDuraMin(int v) { duraMin_ = v; }
	inline void setXOrigin(int v){ xOrigin_ = v; }
	inline void setMinTimeOrigin(int v){ minTimeOrigin_ = v; }
	inline void setMinYTimeOrigin(int v){ minYTimeOrigin_ = v; }
	inline void setMaxYTimeOrigin(int v){ maxYTimeOrigin_ = v; }
	inline void setPersisThresh(double v){ persistenceThreshold_ = v; }
	inline void setMinSeg(std::vector<ttk::SimplexId> &v){ minSeg_ = &v; }
	inline void setSaddleSeg(std::vector<ttk::SimplexId> &v){ saddleSeg_ = &v; }
	inline void setErrSurf(double v){ errSurf_ = v;}
	inline void setOnlyFrameSurface(bool v){ onlyFrameSurface_ = v;}
	inline void setMaxSurfSize(int m){ maxSurfSize_ = m;}
	inline void setBoundaryX(double m){boundaryX_ = m;}
	inline void setBoundaryXMin_(double m){boundaryXMin_ = m;}
	inline void setBoundaryYMin_(double m){boundaryYMin_ = m;}

	inline void setBoundaryY(double m){boundaryY_ = m;}
   
   	struct LinearTrajectory {
      double ax, bx, ay, by;
      int startFrame, endFrame;
      int finalChainId = -1;
      double evalX(int t) const { return ax * t + bx; }
      double evalY(int t) const { return ay * t + by; }

	  std::vector<std::pair<int, ttk::SimplexId>> criticalPoints;

      ttk::SimplexId getOriginalVertex(int frame) const {
        for(const auto &cp : criticalPoints) {
          if(cp.first == frame) return cp.second;
        }
        return -1;
      }	  
    };


    struct FuseRecord {
      int i, j;           
      int endFrame;       
      int startFrame;     
      int finalContrib;
    };


    template <class dataType, class triangulationType>
    int execute(         
                std::vector<int> &durations,           
                std::vector<double> &VX,
                std::vector<double> &VY,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMean,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int frameSurf,
                std::vector<LinearTrajectory> &merge,
                const triangulationType *triangulation);

    
    int correctTrajectory(
        std::vector<std::vector<int>>    &trajTime,
        std::vector<std::vector<int>>    &trajVertexId,
        std::vector<std::vector<double>> &coordsX,
        std::vector<std::vector<double>> &coordsY,
        std::vector<LinearTrajectory> &merge,
        std::vector<LinearTrajectory> &newTraj,
        std::vector<FuseRecord> &fuseRecords
    );

   protected: 



    #ifdef TTK_ENABLE_EIGEN
    int linearRegression(
        const std::vector<int> &T,   
        const std::vector<double> &X,  
        const std::vector<double> &Y, 
        LinearTrajectory &traj
    );
    #endif



    int computeMeanUnitDirectionLinear(
      const std::vector<LinearTrajectory> &newTraj,
      std::vector<double> &meanDx,
      std::vector<double> &meanDy,
      std::vector<double> &meanDz
    );


	template <class dataType, class triangulationType>
	int computeMergeTree(
               const ttk::SimplexId frameSurf,
               const triangulationType *triangulation,
               std::vector<LinearTrajectory> &finalTraj,
               std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
               std::vector<double>              &surfMin,
               std::vector<double>              &surfMax,
               std::vector<double>              &surfMean
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
    double maxRadius_;
    int maxFrameDist_;
    double spatialScale_;
    double interFrame_;
    bool convertDur_;
	bool onlyFrameSurface_;
    double minVx_;
	double maxVx_;
	bool enableFilteringMinVx_;
	bool enableFilteringTimeOrigin_;
	bool enableFilteringCosY_;
	bool enableFilteringDuration_;
	int duraMin_;
	int xOrigin_;
	int minTimeOrigin_;
	int minYTimeOrigin_;
	int maxYTimeOrigin_;
	int minFrameDist_;
	double boundaryY_;
	double boundaryYMin_;
	double boundaryXMin_;
	double boundaryX_;
	double persistenceThreshold_;
	std::vector<ttk::SimplexId> *minSeg_;
	std::vector<ttk::SimplexId> *saddleSeg_;
	double errSurf_;
	int maxSurfSize_;
  }; // DebrisTracer class

} // namespace ttk



#ifdef TTK_ENABLE_EIGEN
int ttk::DebrisTracer::linearRegression(
  const std::vector<int> &T,   // times t_i
  const std::vector<double> &X,   // positions x_i
  const std::vector<double> &Y,   // positions y_i
  LinearTrajectory &traj
) {
  const int n = static_cast<int>(T.size());
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
  traj.ax = bxv[0];
  traj.bx = bxv[1];
  traj.ay = byv[0];
  traj.by = byv[1];
  traj.finalChainId = -1;

  return 1;
}
#endif


/**
 * Count the number of unique VTK cell ids incident to the given
 * vertex set (via vertex-star traversal). 
 */
int ttk::DebrisTracer::computeSurfaceCellCount(
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
int ttk::DebrisTracer::computeMeanUnitDirectionLinear(
  const std::vector<LinearTrajectory> &newTraj,
  std::vector<double> &meanDx,
  std::vector<double> &meanDy,
  std::vector<double> &meanDz
)  {
  const size_t nTraj = newTraj.size();
  meanDx.assign(nTraj, 0.0);
  meanDy.assign(nTraj, 0.0);
  meanDz.assign(nTraj, 0.0);

  for(size_t i = 0; i < nTraj; ++i) {
    const auto &t = newTraj[i];
    double mag = std::sqrt(t.ax *t.ax  + t.ay*t.ay + 1.0);
    if(mag > 0.0) {
      meanDx[i] = t.ax / mag;
      meanDy[i] = t.ay / mag;
      meanDz[i] = 1.0 / mag;
    }
  }
  return 1;
}

/**
 * Fuse and linearize trajectories into longer segments under
 * direction, temporal-gap and spatial-distance constraints.
 */
int ttk::DebrisTracer::correctTrajectory(
    std::vector<std::vector<int>>    &trajTime,
    std::vector<std::vector<int>>    &trajVertexId,
    std::vector<std::vector<double>> &coordsX,
    std::vector<std::vector<double>> &coordsY,
    std::vector<LinearTrajectory> &merge,
    std::vector<LinearTrajectory> &newTraj,
    std::vector<FuseRecord> &fuseRecords
){
  ttk::Timer timer;
  const int numTraj = static_cast<int>(trajTime.size());
  this->printMsg("Linearization and chaining (" + std::to_string(numTraj)
                 + " input trajectories)");

  const double x_min = boundaryXMin_;
  const double x_max = boundaryX_;
  const double y_min = boundaryYMin_;
  const double y_max = boundaryY_;

  auto dirDot = [&](int i, int j,
                    const std::vector<double> &meanDx,
                    const std::vector<double> &meanDy,
                    const std::vector<double> &meanDz) -> double {
    return meanDx[i] * meanDx[j] + meanDy[i] * meanDy[j] + meanDz[i] * meanDz[j];
  };

  auto temporalOk = [&](int sFrame, int eFrame) -> bool {
    return (sFrame - eFrame > minFrameDist_) && (sFrame - eFrame < maxFrameDist_);
  };

  auto dist2AtStartFrame = [&](const LinearTrajectory &coefI,
                               const LinearTrajectory &coefJ,
                               int sFrame) -> double {
    const double xTh = coefI.evalX(sFrame);
    const double yTh = coefI.evalY(sFrame);
    const double xJ  = coefJ.evalX(sFrame);
    const double yJ  = coefJ.evalY(sFrame);
    const double dx = xJ - xTh, dy = yJ - yTh;
    return dx * dx + dy * dy;
  };

  auto resetContribChain = [&](const std::vector<FuseRecord> &chain) {
    for(const auto &fr : chain) {
      newTraj[fr.i].finalChainId = -1;
      newTraj[fr.j].finalChainId = -1;
    }
  };

  auto violatesBBox = [&](const LinearTrajectory &c) -> bool {
    if(!enableFilteringTimeOrigin_) return false;
    const double autoMinX = xOrigin_ - 0.3 * boundaryX_;
    const double autoMaxX = xOrigin_ + 0.3 * boundaryX_;
    return (c.bx < autoMinX || c.bx > autoMaxX || c.by > boundaryY_ || c.by < boundaryXMin_);
  };

  auto passInclinationYNy = [&](double ny) -> bool {
    if(!enableFilteringCosY_) { return true; }
    return (-filtreY_ <= ny && ny <= filtreY_);
  };

  auto passSpeedXAx = [&](double ax) -> bool {
    const double vx_abs = ax * spatialScale_ * (1.0 / interFrame_);
    return (!enableFilteringMinVx_) ? true : (vx_abs >= minVx_ && vx_abs <= maxVx_);
  };

  auto passDirSpeedNyAx = [&](double ny, double ax) -> bool {
    if(!enableFilteringMinVx_ && !enableFilteringCosY_) {
      return true;
    }
    if(!passInclinationYNy(ny)) {
      return false;
    }
    return passSpeedXAx(ax);
  };

  auto passDirSpeed = [&](const LinearTrajectory &c) -> bool {
    const double mag = std::sqrt(c.ax * c.ax + c.ay * c.ay + 1.0);
    const double ny  = c.ay / mag;
    return passDirSpeedNyAx(ny, c.ax);
  };

  auto passDura = [&](const LinearTrajectory &c) -> bool {
    if (!enableFilteringDuration_) { return true; }
    return (duraMin_ <= std::abs(c.endFrame - c.startFrame));
  };

  auto passTimeOrigin = [&](const LinearTrajectory &c) -> bool {
    if (!enableFilteringTimeOrigin_) return true; 
    if(std::abs(c.ax) < 1e-8) return true;
    const double tCross = (xOrigin_ - c.bx) / c.ax;
    const double yCross = c.ay * tCross + c.by;
    return (tCross >= minTimeOrigin_ && yCross >= minYTimeOrigin_ && yCross <= maxYTimeOrigin_);
  };

  auto inFinalBox = [&](double x, double y) -> bool {
    return (x >= x_min && x <= x_max && y >= y_min && y <= y_max);
  };


  auto buildSamplesForChain = [&](const std::vector<FuseRecord> &chain,
                                  std::vector<int> &T, std::vector<double> &X, std::vector<double> &Y) {
    int capacity = static_cast<int>(chain.size()) * 2 + 2;
    T.reserve(capacity); X.reserve(capacity); Y.reserve(capacity);
    for(const auto &r : chain) {
      std::vector<int> T2{trajTime[r.i].front(), r.endFrame};
      const auto &cI = newTraj[r.i];
      for(const int t : T2) {
        X.push_back(cI.evalX(t));
        Y.push_back(cI.evalY(t));
        T.push_back(t);
      }
    }
    const FuseRecord &r = chain.back();
    const auto &cJ = newTraj[r.j];
    X.push_back(cJ.evalX(r.startFrame));
    X.push_back(cJ.evalX(trajTime[r.j].back()));
    Y.push_back(cJ.evalY(r.startFrame));
    Y.push_back(cJ.evalY(trajTime[r.j].back()));
    T.push_back(r.startFrame);
    T.push_back(trajTime[r.j].back());
  };

  auto fitLineCoefForChain = [&](const std::vector<FuseRecord> &chain) -> LinearTrajectory {
    std::vector<int>    T;
    std::vector<double> X, Y;
    buildSamplesForChain(chain, T, X, Y);
    LinearTrajectory lineCoef;
    linearRegression(T, X, Y, lineCoef);
    lineCoef.startFrame = trajTime[chain[0].i].front();
    lineCoef.endFrame   = trajTime[chain.back().j].back();
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
  const double maxLinkDist2        = maxRadius_;

  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || trajTime[i].empty()) continue;

    const int endFrame = trajTime[i].back();

    double bestDot   = similarityThreshold;
    double bestDist2 = std::numeric_limits<double>::infinity();
    int    bestJ     = -1;
	double bestTime = maxFrameDist_; 

    for(int j = 0; j < numTraj; ++j) {
      if(usedAsEnd[j] || j == i || trajTime[j].empty()) continue;

      const int startFrame = trajTime[j].front();

      const double dist2 = dist2AtStartFrame(newTraj[i], newTraj[j], startFrame);
      if(dist2 > maxLinkDist2) continue;

      const double dot = dirDot(i, j, meanDx, meanDy, meanDz);
      if(dot < bestDot) continue;

	  if(!temporalOk(startFrame, endFrame)) continue;
	  if (std::abs(endFrame-startFrame) > bestTime) continue;

      if(dist2 < bestDist2) {
        bestDot   = dot;
        bestDist2 = dist2;
        bestJ     = j;
		bestTime = endFrame - startFrame;
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

    std::vector<FuseRecord> chain{r1};   
   	r1.finalContrib = finalId;
    newTraj[r1.i].finalChainId = finalId;
    newTraj[r1.j].finalChainId = finalId;
    used[idx1] = true;

    // prepend
    bool prepended = true;
    while(prepended) {
      prepended = false;
      for(size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
        if(used[idx2]) continue;
        auto &r2 = fuseRecords[idx2];
        if(r2.j == chain.front().i) {
          chain.insert(chain.begin(), r2);
          r2.finalContrib = finalId;
          newTraj[r2.i].finalChainId = finalId;
          newTraj[r2.j].finalChainId = finalId;
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
        if(chain.back().j == r2.i) {
          chain.push_back(r2);
          r2.finalContrib = finalId;
          newTraj[r2.i].finalChainId = finalId;
          newTraj[r2.j].finalChainId = finalId;
          used[idx2] = true;
          extended = true;
          break;
        }
      }
    }

    LinearTrajectory lineCoef = fitLineCoefForChain(chain);

    if(passDirSpeed(lineCoef) && passTimeOrigin(lineCoef) && passDura(lineCoef)) {
      if(violatesBBox(lineCoef)) {
        resetContribChain(chain);
      } else {
        { // add initials criticalPoints 
          const int firstTraj = chain[0].i;
          for(size_t k = 0; k < trajTime[firstTraj].size(); ++k) {
            lineCoef.criticalPoints.emplace_back(
              trajTime[firstTraj][k],
              static_cast<ttk::SimplexId>(trajVertexId[firstTraj][k])
            );
          }
          for(const auto &rec : chain) {
            const int tj = rec.j;
            for(size_t k = 0; k < trajTime[tj].size(); ++k) {
              lineCoef.criticalPoints.emplace_back(
                trajTime[tj][k],
                static_cast<ttk::SimplexId>(trajVertexId[tj][k])
              );
            }
          }
		}
        merge.push_back(lineCoef);
      }
    } else {
      resetContribChain(chain);
    }
  }


  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || usedAsEnd[i] || trajTime[i].empty()) continue;

    const double nyOrphan = meanDy[i];
    const double axOrphan = newTraj[i].ax;
    if(passDirSpeedNyAx(nyOrphan, axOrphan)) {

      LinearTrajectory lineCoef = newTraj[i];
      lineCoef.startFrame = trajTime[i].front();
      lineCoef.endFrame   = trajTime[i].back();

      if(violatesBBox(lineCoef)) {
        continue;
      }
      if(!passTimeOrigin(lineCoef) || !passDura(lineCoef)) { continue; }
      lineCoef.criticalPoints.reserve(trajTime[i].size());
      for(size_t k = 0; k < trajTime[i].size(); ++k) {
        lineCoef.criticalPoints.emplace_back(
          trajTime[i][k],
          static_cast<ttk::SimplexId>(trajVertexId[i][k])
        );
      }
      newTraj[i].finalChainId = static_cast<int>(merge.size());
      merge.push_back(lineCoef);
    }
  }
  

  for(auto &c : merge) {
    int start = c.startFrame;
    int end   = c.endFrame;

    if(end < start) {
      std::swap(c.startFrame, c.endFrame);
      end = c.endFrame;
      start = c.startFrame;
    }

    while(end > start) {
      const double xEnd = c.evalX(end);
      const double yEnd = c.evalY(end);

      if(inFinalBox(xEnd, yEnd)) {
        break;
      }
      --end;
    }

    c.endFrame = end;
  }

  this->printMsg("Linearization and chaining (" + std::to_string(merge.size())
                 + " output chains)", 1.0, timer.getElapsedTime(),
                 this->threadNumber_);
  return 1;
}


template <class dataType, class triangulationType>
int ttk::DebrisTracer::execute(
                std::vector<int>                &durations,
                std::vector<double>             &VX,
                std::vector<double>             &VY,
                std::vector<double>             &surfMin,
                std::vector<double>             &surfMax,
                std::vector<double>             &surfMean,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int frameSurf,
                std::vector<LinearTrajectory> &finalTraj,
                const triangulationType *triangulation) {

    ttk::Timer timer;
    this->printMsg("Computing statistics for "
                   + std::to_string(finalTraj.size()) + " trajectories");
    const int numMerge = finalTraj.size();

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    
    for(int i = 0; i < numMerge; ++i) {
        if (convertDur_) durations[i] = (finalTraj[i].endFrame - finalTraj[i].startFrame)*interFrame_;
        else durations[i] = finalTraj[i].endFrame - finalTraj[i].startFrame;
    } 
    
	double conversion = spatialScale_*(1/interFrame_);
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int i = 0; i < numMerge; ++i) {
        VX[i] = finalTraj[i].ax * conversion;
        VY[i] = finalTraj[i].ay * conversion;
    }
	computeMergeTree<dataType, triangulationType>(
			frameSurf,
			triangulation,
			finalTraj,
			allVertexDebris,
			surfMin, surfMax, surfMean);


    this->printMsg("Statistics complete", 1.0, timer.getElapsedTime(),
                   this->threadNumber_);
    return 1;
}

            

template <class dataType, class triangulationType>
int ttk::DebrisTracer::computeMergeTree(
  const ttk::SimplexId frameSurf,
  const triangulationType *triangulation,
  std::vector<LinearTrajectory> &finalTraj,
  std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
  std::vector<double>              &surfMin,
  std::vector<double>              &surfMax,
  std::vector<double>              &surfMean
) {

  ttk::Timer globalTimer;
  const ttk::SimplexId nPixels = triangulation->getNumberOfVertices();
  const int nFrames = (!onlyFrameSurface_) ? inputData_.size() : 1;
  std::vector<std::vector<double>> trajSurfaces(finalTraj.size());
  const auto nTraj = finalTraj.size();
  this->printMsg("Merge-tree surface segmentation (" + std::to_string(nFrames)
                 + " frame(s), " + std::to_string(nPixels) + " vertices)");
  std::vector<char> trajDouble(nTraj);
  for(int frame = 0  ; frame < nFrames; frame++) {

    ttk::Timer frameTimer;
    std::fill(trajDouble.begin(), trajDouble.end(), 0);
  	frame = (!onlyFrameSurface_)  ? frame : frameSurf;
	this->printMsg("Processing frame " + std::to_string(frame) + "/"
                   + std::to_string(nFrames - 1));
    auto *scalars = static_cast<dataType *>(inputData_[frame]);

    std::vector<ttk::SimplexId> pdOffsets(nPixels);
    ttk::preconditionOrderArray<dataType>(
      static_cast<size_t>(nPixels),
      scalars,
      pdOffsets.data(),
      this->threadNumber_);

    ttk::PersistenceDiagram persistenceDiagram;
    persistenceDiagram.setThreadNumber(this->threadNumber_);
    persistenceDiagram.setDebugLevel(0);
    persistenceDiagram.setBackend(
      ttk::PersistenceDiagram::BACKEND::DISCRETE_MORSE_SANDWICH);
    persistenceDiagram.preconditionTriangulation(
      const_cast<triangulationType *>(triangulation));

    ttk::DiagramType diagram;

    const int statusPD = persistenceDiagram.execute(
      diagram,
      scalars,
      0, 
      pdOffsets.data(),
      const_cast<triangulationType *>(triangulation));

    if(statusPD != 0) {
      this->printErr("PersistenceDiagram::execute failed");
      return -1;
    }

    // Critical Points tresh 
	double maxPers = 0.0;

	for(const auto &pair : diagram) {
	  if(pair.dim != 0) continue;
	  const double pers = pair.persistence();
	  if(!std::isfinite(pers)) continue;
	  if(pers > maxPers) maxPers = pers;
	}

	const double threshold = maxPers * (this->persistenceThreshold_ / 100.0);

	std::vector<ttk::SimplexId> criticalPoints;
	criticalPoints.reserve(diagram.size() * 2);

	for(const auto &pair : diagram) {
	  if(pair.dim != 0) continue;
	  const double pers = pair.persistence();
	  if(!std::isfinite(pers)) continue;
	  if(pers < threshold) continue;
	  criticalPoints.push_back(pair.birth.id);
	  criticalPoints.push_back(pair.death.id);
	}

    // Topological Simplification
    std::vector<dataType> outScalars(nPixels);
    std::copy(scalars, scalars + nPixels, outScalars.begin());

    std::vector<ttk::SimplexId> offsets = pdOffsets;

    ttk::TopologicalSimplification topoSimp;
    topoSimp.setThreadNumber(this->threadNumber_);
    topoSimp.setDebugLevel(0);
    topoSimp.setBackend(ttk::TopologicalSimplification::BACKEND::LTS);
    topoSimp.preconditionTriangulation(const_cast<triangulationType *>(triangulation));

    const bool addPerturbation = true;
    const ttk::SimplexId constraintNumber = static_cast<ttk::SimplexId>(criticalPoints.size());
    const ttk::DiagramType emptyDiagram;
    topoSimp.execute<dataType, triangulationType>(
      scalars,
      outScalars.data(),
      criticalPoints.empty() ? nullptr : criticalPoints.data(),
      pdOffsets.data(),
      offsets.data(),
      constraintNumber,
      addPerturbation,
      *const_cast<triangulationType *>(triangulation),
      emptyDiagram);

    // Merge tree (ttk::ExTreeM)

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
    pathComp.setDebugLevel(0);
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
    exTreeM.setDebugLevel(0);

    std::vector<ttk::SimplexId> orderJoin(order);

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

    //  Association to trajectories
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

    std::vector<char> segCleaned(segmentId.size(), 0);

    for(size_t trajId = 0; trajId < nTraj; ++trajId) {
      const auto &traj = finalTraj[trajId];

      if(frame < traj.startFrame || frame > traj.endFrame)
        continue;

      ttk::SimplexId vId = -1;
	  vId = traj.getOriginalVertex(frame);
      if(vId<0) {
        const double x = traj.evalX(frame);
        if(x < boundaryXMin_ || x > boundaryX_+1)
          continue;

        const double y = traj.evalY(frame);
        if(y < boundaryYMin_ || y > boundaryY_+1)
          continue;

        const ttk::SimplexId xi = std::lround(x);
        const ttk::SimplexId yi = std::lround(y);
        vId = xi + yi * (boundaryX_ - boundaryXMin_ + 1);
      }

      if(vId < 0 || vId >= nPixels)
        continue;

	  if(regionType[vId] == 0) {

	    auto segId = segmentation[vId];

		if(segId >= 0 && segId < (ttk::SimplexId)segmentId.size() && !segCleaned[segId] && segmentId[segId].size() > 8 && errSurf_ != 0) 		 {
	  	  cleanDarkSegmentInPlace<dataType, triangulationType>(
	  	    segmentId[segId], scalars, triangulation, static_cast<int>(errSurf_));
	  	  segCleaned[segId] = 1;
	    }

		if (static_cast<int>(segmentId[segId].size()) > maxSurfSize_) continue;

	    double surfVal = static_cast<double>(
	  	computeSurfaceCellCount(segmentId[segId], triangulation));
	    if(surfVal == 0) surfVal = 1;

	    trajSurfaces[trajId].push_back(surfVal);
		
		for (size_t i = 0; i<segmentId[segId].size(); i++){
			int v = segmentId[segId][i];
			int check = allVertexDebris[frame][v];
		  	if (check == -1 ) allVertexDebris[frame][v] = trajId;
			else if (check != static_cast<int>(trajId) ){
				trajDouble[trajId] = 1;
				if (check >=0)
					trajDouble[check]=1;
			}
		}

	    (*saddleSeg_)[trajId]   = segMinVertex[segId];
	  	(*minSeg_)[trajId]      = segId;
	  }

    }
    for(int t = 0; t < static_cast<int>(trajDouble.size()); t++) {
      if(!trajDouble[t]) continue;
      for(int v = 0; v < nPixels; v++) {
        if(allVertexDebris[frame][v] == t)
          allVertexDebris[frame][v] = -2;
      }
    }
    this->printMsg("Frame " + std::to_string(frame) + " done", 1.0,
                   frameTimer.getElapsedTime(), this->threadNumber_);
  }

  // Per-trajectory surface statistics (surfMin, surfMax, surfMean)
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
      surfMean[trajId]   = mean;
    } else {
      surfMin[trajId] = 0.0;
      surfMax[trajId] = 0.0;
      surfMean[trajId] = 0.0;
    }
  }

  this->printMsg("Merge-tree surface segmentation complete", 1.0,
                 globalTimer.getElapsedTime(), this->threadNumber_);
  return 0;
}


template <class dataType>
dataType ttk::DebrisTracer::otsuThresholdLocal(
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

template <typename dataType, typename triangulationType>
void ttk::DebrisTracer::cleanDarkSegmentInPlace(
								   std::vector<ttk::SimplexId> &segmentVerts,
                                   const dataType *scalars,
                                   const triangulationType *triangulation,
                                   const int otsuBins) {
  if(segmentVerts.size() < 2)
    return;

  const dataType T = otsuThresholdLocal<dataType>(segmentVerts, scalars, otsuBins);
  std::vector<char> inSeg(triangulation->getNumberOfVertices(), 0);
  std::vector<ttk::SimplexId> kept;
  kept.reserve(segmentVerts.size());

  for(const auto v : segmentVerts) {
    if(scalars[v] <= T) {
      kept.push_back(v);
      inSeg[v] = 1;
    }
  }

  const size_t segSize = segmentVerts.size();
  const size_t minKeep = std::max<size_t>(2, (size_t)std::ceil(0.10 * (double)segSize));   if(kept.size() < minKeep && segSize >= 2) {
  
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
