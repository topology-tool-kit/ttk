/// \ingroup base
/// \class ttk::TrackingPostProcessing
/// \author Théophane Loloum <theophane.loloum@gmail.com>
/// \date April 2026
///
/// \brief TTK processing package for post-processing tracked trajectories:
/// linearization, fusion (chaining), and merge-tree-based segmentation
/// statistics.
///
/// This module takes the initial trajectories produced by an upstream tracker
/// (TrackingFromFields) and refines them:
///   - linearization: fit a 2D line (x(t)=ax*t+bx, y(t)=ay*t+by) through each
///     trajectory point cloud via least squares (Eigen);
///   - fusion: greedily chain temporally-adjacent, directionally-consistent
///     linearized segments into longer trajectories;
///   - merge-tree segmentation: per-frame, compute a merge tree of the scalar
///     field and associate each trajectory
///
/// \sa ttk::TrackingFromFields
/// \sa ttk::TrackingFromCriticalPoints

#pragma once

#include <DataTypes.h>
#include <Debug.h>
#include <Geometry.h>
#include <Timer.h>
#include <Triangulation.h>

// merge-tree
#include <ExTreeM.h>
#include <FTMTreePP.h>
#include <LocalizedTopologicalSimplification.h>
#include <PathCompression.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ttk {

  class TrackingPostProcessing : virtual public Debug {

  public:
    /// @brief Linear trajectory: x(t) = ax*t + bx, y(t) = ay*t + by,
    /// defined on the inclusive frame range [startFrame, endFrame].
    struct LinearTrajectory {
      double ax{0.0}, bx{0.0}, ay{0.0}, by{0.0};
      int startFrame{0};
      int endFrame{0};
      int finalChainId{-1};
      int originalTrajId{-1};
      bool isLinearized{true};

      std::vector<std::pair<int, ttk::SimplexId>> criticalPoints;

      inline double evalX(double t) const {
        return ax * t + bx;
      }
      inline double evalY(double t) const {
        return ay * t + by;
      }

      inline ttk::SimplexId getOriginalVertex(int frame) const {
        for(const auto &cp : criticalPoints) {
          if(cp.first == frame)
            return cp.second;
        }
        return -1;
      }
    };

    /// @brief Fusion-link record: trajectory i ends and trajectory j starts
    struct FuseRecord {
      int i{-1}, j{-1};
      int endFrame{0};
      int startFrame{0};
      int finalContrib{-1};
    };

    TrackingPostProcessing();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      return triangulation->preconditionVertexStars();
    }

    inline void setInputScalars(const std::vector<void *> &inputScalars) {
      inputData_ = inputScalars;
    }
    inline void setCosCol(double v) {
      cosCol_ = v;
    }
    inline void setMaxRadius(double v) {
      maxRadius_ = v;
    }
    inline void setMaxFrameDist(int v) {
      maxFrameDist_ = v;
      minFrameDist_ = -v;
    }

    inline void setPersistenceThreshold(double v) {
      persistenceThreshold_ = v;
    }
    inline void setMaxSurfSize(int v) {
      maxSurfSize_ = v;
    }
    inline void setUseOtsuSimplification(bool v) {
      useOtsuSimplification_ = v;
    }
    inline void setOtsuBins(int v) {
      otsuBins_ = v;
    }

    inline void setBoundaryXMin(double v) {
      boundaryXMin_ = v;
    }
    inline void setBoundaryXMax(double v) {
      boundaryXMax_ = v;
    }
    inline void setBoundaryYMin(double v) {
      boundaryYMin_ = v;
    }
    inline void setBoundaryYMax(double v) {
      boundaryYMax_ = v;
    }

    inline void setDoLinearize(bool v) {
      doLinearize_ = v;
    }
    inline void setDoFusion(bool v) {
      doFusion_ = v;
    }
    inline void setDoLinearizeFuse(bool v) {
      doLinearizeFuse_ = v;
    }
    inline void setDoMergeTree(bool v) {
      doMergeTree_ = v;
    }
    inline void setUseSplitTree(int v) {
      useSplitTree_ = v;
    }

    /// @brief Linearize + (optional) chain input per-trajectory point clouds.
    ///
    /// @param[in]  trajTime      per-trajectory, frame indices (sorted)
    /// @param[in]  trajVertexId  per-trajectory, vertex global ids
    /// @param[in]  coordsX       per-trajectory, X coordinates
    /// @param[in]  coordsY       per-trajectory, Y coordinates
    /// @param[out] linearTraj    per-input-trajectory linear fits (with
    ///                           finalChainId set when the segment is part of
    ///                           a fused chain in outputTraj, -1 otherwise)
    /// @param[out] outputTraj    final trajectory set (fused chains + un-fused
    ///                           survivors)
    /// @param[out] fuseRecords   list of i -> j fusion links
    int correctTrajectory(const std::vector<std::vector<int>> &trajTime,
                          const std::vector<std::vector<int>> &trajVertexId,
                          const std::vector<std::vector<double>> &coordsX,
                          const std::vector<std::vector<double>> &coordsY,
                          const std::vector<int> &trajCriticalType,
                          std::vector<LinearTrajectory> &linearTraj,
                          std::vector<LinearTrajectory> &outputTraj,
                          std::vector<FuseRecord> &fuseRecords);

    /// @brief Compute merge-tree-based segmentation per trajectory && per
    /// frame.
    ///
    /// @param[in]  triangulation        triangulation of the scalar field
    /// @param[in,out] finalTraj         trajectories to annotate
    /// @param[out] surfMin              per-trajectory minimum surface
    /// @param[out] surfMax              per-trajectory maximum surface
    /// @param[out] surfMean             per-trajectory mean surface
    /// @param[out] vertexTrajPerFrame   per-frame, per-vertex labelling:
    ///                                  -1 = no surface, -2 = collision,
    template <class dataType, class triangulationType>
    int computeMergeTree(const triangulationType *triangulation,
                         const std::vector<LinearTrajectory> &finalTraj,
                         std::vector<double> &surfMin,
                         std::vector<double> &surfMax,
                         std::vector<double> &surfMean,
                         std::vector<std::vector<int>> &vertexTrajPerFrame);

    /// @brief correctTrajectory + computeMergeTree.
    template <class dataType, class triangulationType>
    int execute(const std::vector<std::vector<int>> &trajTime,
                const std::vector<std::vector<int>> &trajVertexId,
                const std::vector<std::vector<double>> &coordsX,
                const std::vector<std::vector<double>> &coordsY,
                const std::vector<int> &trajCriticalType,
                std::vector<LinearTrajectory> &linearTraj,
                std::vector<LinearTrajectory> &finalTraj,
                std::vector<FuseRecord> &fuseRecords,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMean,
                std::vector<std::vector<int>> &vertexTrajPerFrame,
                const triangulationType *triangulation);

  protected:
#ifdef TTK_ENABLE_EIGEN
    int linearRegression(const std::vector<int> &T,
                         const std::vector<double> &X,
                         const std::vector<double> &Y,
                         LinearTrajectory &traj);
#endif

    int computeMeanUnitDirectionLinear(
      const std::vector<LinearTrajectory> &newTraj,
      std::vector<std::array<double, 3>> &meanDir);

    int
      computeSurfaceCellCount(const std::vector<ttk::SimplexId> &surfVertices,
                              const ttk::AbstractTriangulation *triangulation);

    template <class dataType>
    dataType otsuThresholdLocal(const std::vector<ttk::SimplexId> &verts,
                                const dataType *scalars,
                                const int nbins);

    template <class dataType, class triangulationType>
    void cleanDarkSegmentInPlace(std::vector<ttk::SimplexId> &segmentVerts,
                                 const dataType *scalars,
                                 const triangulationType *triangulation,
                                 const int otsuBins);

    std::vector<void *> inputData_{};

    double cosCol_{0.9};
    double maxRadius_{225.0}; // squared pixel distance
    int maxFrameDist_{30};
    int minFrameDist_{-30};

    double persistenceThreshold_{0.0};
    int maxSurfSize_{10000};
    bool useOtsuSimplification_{false};
    int otsuBins_{0};

    double boundaryXMin_{0.0};
    double boundaryXMax_{0.0};
    double boundaryYMin_{0.0};
    double boundaryYMax_{0.0};

    bool doLinearize_{true};
    bool doFusion_{true};
    bool doLinearizeFuse_{true};
    bool doMergeTree_{false};
    int useSplitTree_{2};
  };

} // namespace ttk

#ifdef TTK_ENABLE_EIGEN
inline int
  ttk::TrackingPostProcessing::linearRegression(const std::vector<int> &T,
                                                const std::vector<double> &X,
                                                const std::vector<double> &Y,
                                                LinearTrajectory &traj) {
  const int n = static_cast<int>(T.size());
  if(n < 1)
    return 0;
  Eigen::MatrixXd M(n, 2);
  Eigen::VectorXd vx(n), vy(n);
  for(int i = 0; i < n; ++i) {
    M(i, 0) = T[i];
    M(i, 1) = 1.0;
    vx(i) = X[i];
    vy(i) = Y[i];
  }
  const Eigen::Vector2d bxv
    = (M.transpose() * M).ldlt().solve(M.transpose() * vx);
  const Eigen::Vector2d byv
    = (M.transpose() * M).ldlt().solve(M.transpose() * vy);
  traj.ax = bxv[0];
  traj.bx = bxv[1];
  traj.ay = byv[0];
  traj.by = byv[1];
  traj.finalChainId = -1;
  return 1;
}
#endif

inline int ttk::TrackingPostProcessing::computeMeanUnitDirectionLinear(
  const std::vector<LinearTrajectory> &newTraj,
  std::vector<std::array<double, 3>> &meanDir) {
  const size_t nTraj = newTraj.size();
  meanDir.assign(nTraj, {0.0, 0.0, 0.0});
  for(size_t i = 0; i < nTraj; ++i) {
    const auto &t = newTraj[i];
    std::array<double, 3> v{t.ax, t.ay, 1.0};
    const double mag = ttk::Geometry::magnitude<double>(v.data(), 3);
    if(mag > 0.0) {
      ttk::Geometry::scaleVector<double>(
        v.data(), 1.0 / mag, meanDir[i].data(), 3);
    }
  }
  return 1;
}

inline int ttk::TrackingPostProcessing::computeSurfaceCellCount(
  const std::vector<ttk::SimplexId> &surfVertices,
  const ttk::AbstractTriangulation *triangulation) {
  std::unordered_set<ttk::SimplexId> cellIds;
  for(const ttk::SimplexId v : surfVertices) {
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

template <class dataType>
dataType ttk::TrackingPostProcessing::otsuThresholdLocal(
  const std::vector<ttk::SimplexId> &verts,
  const dataType *scalars,
  const int nbins) {
  if(verts.empty() || nbins < 2)
    return dataType{0};

  dataType vmin = scalars[verts[0]];
  dataType vmax = scalars[verts[0]];
  for(const auto v : verts) {
    const auto s = scalars[v];
    if(s < vmin)
      vmin = s;
    if(s > vmax)
      vmax = s;
  }
  if(vmax <= vmin)
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
  for(auto &h : hist)
    h /= total;

  std::vector<double> omega(nbins, 0.0);
  std::vector<double> mu(nbins, 0.0);
  omega[0] = hist[0];
  mu[0] = 0.0;
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
  const double t
    = minD + (static_cast<double>(bestK) / (nbins - 1)) * (maxD - minD);
  return static_cast<dataType>(t);
}

template <class dataType, class triangulationType>
void ttk::TrackingPostProcessing::cleanDarkSegmentInPlace(
  std::vector<ttk::SimplexId> &segmentVerts,
  const dataType *scalars,
  const triangulationType *triangulation,
  const int otsuBins) {
  if(segmentVerts.size() < 2)
    return;

  const dataType T
    = otsuThresholdLocal<dataType>(segmentVerts, scalars, otsuBins);
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
  const size_t minKeep
    = std::max<size_t>(2, (size_t)std::ceil(0.10 * (double)segSize));
  if(kept.size() < minKeep && segSize >= 2) {
    for(const auto v : kept)
      inSeg[v] = 0;
    kept.clear();
    std::vector<ttk::SimplexId> tmp = segmentVerts;
    std::sort(tmp.begin(), tmp.end(), [&](ttk::SimplexId a, ttk::SimplexId b) {
      return scalars[a] < scalars[b];
    });
    for(size_t i = 0; i < std::min(minKeep, tmp.size()); ++i) {
      kept.push_back(tmp[i]);
      inSeg[tmp[i]] = 1;
    }
  }
  if(kept.empty())
    return;

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
        if(nb < 0 || !inSeg[nb] || visited[nb])
          continue;
        visited[nb] = 1;
        q.push(nb);
      }
    }
    if(cc.size() > bestCC.size())
      bestCC.swap(cc);
  }
  segmentVerts.swap(bestCC);
}

template <class dataType, class triangulationType>
int ttk::TrackingPostProcessing::computeMergeTree(
  const triangulationType *triangulation,
  const std::vector<LinearTrajectory> &finalTraj,
  std::vector<double> &surfMin,
  std::vector<double> &surfMax,
  std::vector<double> &surfMean,
  std::vector<std::vector<int>> &vertexTrajPerFrame) {

  const size_t nTraj = finalTraj.size();
  surfMin.assign(nTraj, 0.0);
  surfMax.assign(nTraj, 0.0);
  surfMean.assign(nTraj, 0.0);

  if(!doMergeTree_ || nTraj == 0 || inputData_.empty()) {
    vertexTrajPerFrame.clear();
    return 0;
  }

  ttk::Timer globalTimer;
  const ttk::SimplexId nPixels = triangulation->getNumberOfVertices();
  const int nFrames = static_cast<int>(inputData_.size());

  this->printMsg("Merge-tree (" + std::to_string(nFrames) + " f., "
                 + std::to_string(nPixels) + " v., " + std::to_string(nTraj)
                 + " t.)");

  // Per-frame / per-traj accumulated surface, and per-frame collision flags
  std::vector<std::vector<double>> trajSurfPerFrame(
    nFrames, std::vector<double>(nTraj, 0.0));
  std::vector<std::vector<char>> trajDoublePerFrame(
    nFrames, std::vector<char>(nTraj, 0));

  vertexTrajPerFrame.assign(nFrames, std::vector<int>(nPixels, -1));

  int globalError = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(int frame = 0; frame < nFrames; ++frame) {
    if(globalError != 0)
      continue;

    ttk::Timer frameTimer;
    auto *scalars = static_cast<dataType *>(inputData_[frame]);

    // --- Topological Simplification ---
    std::vector<dataType> outScalars(nPixels);
    std::copy(scalars, scalars + nPixels, outScalars.begin());

    std::vector<ttk::SimplexId> offsets(nPixels);
    ttk::preconditionOrderArray<dataType>(
      static_cast<size_t>(nPixels), scalars, offsets.data(), 1);

    dataType sMin = scalars[0], sMax = scalars[0];
    for(ttk::SimplexId i = 1; i < nPixels; ++i) {
      if(scalars[i] < sMin)
        sMin = scalars[i];
      if(scalars[i] > sMax)
        sMax = scalars[i];
    }
    const dataType persThresh = static_cast<dataType>(
      (sMax - sMin) * (this->persistenceThreshold_ / 100.0));

    ttk::lts::LocalizedTopologicalSimplification lts;
    lts.setThreadNumber(1);
    lts.setDebugLevel(0);
    lts.preconditionTriangulation(
      const_cast<triangulationType *>(triangulation));

    const int statusSimp = lts.removeNonPersistentExtrema<
      dataType, ttk::SimplexId, triangulationType>(
      outScalars.data(), offsets.data(), triangulation, persThresh, true,
      ttk::lts::LocalizedTopologicalSimplification::PAIR_TYPE::EXTREMUM_SADDLE);
    if(statusSimp != 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
      globalError = -1;
      continue;
    }

    // --- Merge tree ---
    std::vector<ttk::SimplexId> order(nPixels);
    ttk::preconditionOrderArray<dataType>(
      static_cast<size_t>(nPixels), outScalars.data(), order.data(), 1);

    std::vector<ttk::SimplexId> ascendingManifold(nPixels, -1);
    std::vector<ttk::SimplexId> descendingManifold(nPixels, -1);

    ttk::PathCompression pathComp;
    pathComp.setThreadNumber(1);
    pathComp.setDebugLevel(0);
    pathComp.setComputeSegmentation(true, true, false);

    ttk::PathCompression::OutputSegmentation om{
      ascendingManifold.data(), descendingManifold.data(), nullptr};

    {
      const int statusPC = pathComp.execute<triangulationType>(
        om, order.data(), *const_cast<triangulationType *>(triangulation));
      if(statusPC != 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
        globalError = -1;
        continue;
      }
    }

    std::vector<ttk::SimplexId> orderInv(nPixels);
    for(ttk::SimplexId i = 0; i < nPixels; ++i)
      orderInv[i] = nPixels - order[i] - 1;

    auto &localTrajDouble = trajDoublePerFrame[frame];
    auto &localVertexLabel = vertexTrajPerFrame[frame];
    std::vector<int> vertexTraj(nPixels, -1);

    auto runTree
      = [&](const ttk::SimplexId *mtOrder, ttk::SimplexId *mtManifold,
            ttk::SimplexId *mtScratch) -> bool {
      std::vector<ttk::SimplexId> segmentation(nPixels, -1);
      std::vector<char> regionType(nPixels, 0);
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> persistencePairs;
      std::map<ttk::SimplexId, int> cpMap;
      std::vector<ttk::ExTreeM::Branch> branches;

      ttk::ExTreeM exTreeM;
      exTreeM.setThreadNumber(1);
      exTreeM.setDebugLevel(0);

      const int statusMT = exTreeM.computePairs<triangulationType>(
        persistencePairs, cpMap, branches, segmentation.data(),
        regionType.data(), mtManifold, mtScratch, mtOrder, triangulation,
        ttk::ftm::TreeType::Join);
      if(statusMT != 1)
        return false;

      const ttk::SimplexId maxSegId
        = *std::max_element(segmentation.begin(), segmentation.end());
      std::vector<std::vector<ttk::SimplexId>> segmentId(maxSegId + 1);
      for(size_t vId = 0; vId < segmentation.size(); ++vId) {
        if(regionType[vId] == 0)
          segmentId[segmentation[vId]].push_back(
            static_cast<ttk::SimplexId>(vId));
      }

      std::vector<char> segCleaned(segmentId.size(), 0);

      for(size_t trajId = 0; trajId < nTraj; ++trajId) {
        const auto &traj = finalTraj[trajId];
        if(frame < traj.startFrame || frame > traj.endFrame)
          continue;

        ttk::SimplexId vId = traj.getOriginalVertex(frame);
        if(vId < 0) {
          if(!traj.isLinearized)
            continue;
          const double x = traj.evalX(frame);
          if(x < boundaryXMin_ || x > boundaryXMax_ + 1)
            continue;
          const double y = traj.evalY(frame);
          if(y < boundaryYMin_ || y > boundaryYMax_ + 1)
            continue;
          const ttk::SimplexId xi = static_cast<ttk::SimplexId>(std::lround(x));
          const ttk::SimplexId yi = static_cast<ttk::SimplexId>(std::lround(y));
          vId = xi + yi * (ttk::SimplexId)(boundaryXMax_ - boundaryXMin_ + 1);
        }
        if(vId < 0 || vId >= nPixels)
          continue;
        if(regionType[vId] != 0)
          continue;

        const auto segId = segmentation[vId];
        if(segId < 0 || segId >= (ttk::SimplexId)segmentId.size())
          continue;

        if(trajSurfPerFrame[frame][trajId] > 0.0)
          continue;

        if(useOtsuSimplification_ && !segCleaned[segId]
           && segmentId[segId].size() > 8 && otsuBins_ > 0) {
          cleanDarkSegmentInPlace<dataType, triangulationType>(
            segmentId[segId], scalars, triangulation, otsuBins_);
          segCleaned[segId] = 1;
        }

        if(static_cast<int>(segmentId[segId].size()) > maxSurfSize_)
          continue;

        double surfVal = static_cast<double>(
          computeSurfaceCellCount(segmentId[segId], triangulation));
        if(surfVal == 0)
          surfVal = 1;
        trajSurfPerFrame[frame][trajId] = surfVal;

        const int currentChainId = finalTraj[trajId].finalChainId;

        for(const auto v : segmentId[segId]) {
          const int check = vertexTraj[v];
          if(check == -1) {
            vertexTraj[v] = static_cast<int>(trajId);
          } else if(check != static_cast<int>(trajId)) {
            localTrajDouble[trajId] = 1;
            if(check >= 0)
              localTrajDouble[check] = 1;
          }
        }

        if(currentChainId >= 0) {
          std::unordered_set<ttk::SimplexId> dilatedSet;
          dilatedSet.reserve(segmentId[segId].size() * 4);
          for(const ttk::SimplexId v : segmentId[segId]) {
            const ttk::SimplexId starCount
              = triangulation->getVertexStarNumber(v);
            for(ttk::SimplexId k = 0; k < starCount; ++k) {
              ttk::SimplexId cellId;
              triangulation->getVertexStar(v, k, cellId);
              const int nCellVerts = triangulation->getCellVertexNumber(cellId);
              for(int cv = 0; cv < nCellVerts; ++cv) {
                ttk::SimplexId vDil;
                triangulation->getCellVertex(cellId, cv, vDil);
                dilatedSet.insert(vDil);
              }
            }
          }

          for(const ttk::SimplexId v : dilatedSet) {
            const int prev = localVertexLabel[v];
            if(prev == -1) {
              localVertexLabel[v] = currentChainId;
            } else if(prev != currentChainId && prev != -2) {
              localVertexLabel[v] = -2;
            }
          }
        }
      } // trajectory loop
      return true;
    };

    bool ok = true;
    if(useSplitTree_ == 1) {
      ok = runTree(
        order.data(), descendingManifold.data(), ascendingManifold.data());
    } else if(useSplitTree_ == 0) {
      ok = runTree(
        orderInv.data(), ascendingManifold.data(), descendingManifold.data());
    } else {
      ok = runTree(
        order.data(), descendingManifold.data(), ascendingManifold.data());
      if(ok)
        ok = runTree(
          orderInv.data(), ascendingManifold.data(), descendingManifold.data());
    }

    if(!ok) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
      globalError = -1;
      continue;
    }
  } // frame loop

  if(globalError != 0) {
    this->printErr("Error in merge-tree frame processing");
    return -1;
  }

  for(int frame = 0; frame < nFrames; ++frame) {
    for(size_t trajId = 0; trajId < nTraj; ++trajId) {
      if(trajDoublePerFrame[frame][trajId])
        trajSurfPerFrame[frame][trajId] = 0.0;
    }
  }

  // Per-trajectory stat
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t trajId = 0; trajId < nTraj; ++trajId) {
    double minVal = std::numeric_limits<double>::max();
    double maxVal = 0.0;
    double sum = 0.0;
    int count = 0;
    for(int frame = 0; frame < nFrames; ++frame) {
      const double s = trajSurfPerFrame[frame][trajId];
      if(s > 0.0) {
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
      surfMean[trajId] = sum / static_cast<double>(count);
    }
  }

  this->printMsg("Segmentation complete", 1.0, globalTimer.getElapsedTime(),
                 this->threadNumber_);
  return 0;
}

template <class dataType, class triangulationType>
int ttk::TrackingPostProcessing::execute(
  const std::vector<std::vector<int>> &trajTime,
  const std::vector<std::vector<int>> &trajVertexId,
  const std::vector<std::vector<double>> &coordsX,
  const std::vector<std::vector<double>> &coordsY,
  const std::vector<int> &trajCriticalType,
  std::vector<LinearTrajectory> &linearTraj,
  std::vector<LinearTrajectory> &finalTraj,
  std::vector<FuseRecord> &fuseRecords,
  std::vector<double> &surfMin,
  std::vector<double> &surfMax,
  std::vector<double> &surfMean,
  std::vector<std::vector<int>> &vertexTrajPerFrame,
  const triangulationType *triangulation) {

  ttk::Timer timer;

  this->correctTrajectory(trajTime, trajVertexId, coordsX, coordsY,
                          trajCriticalType, linearTraj, finalTraj, fuseRecords);

  const int numFinal = static_cast<int>(finalTraj.size());
  surfMin.assign(numFinal, 0.0);
  surfMax.assign(numFinal, 0.0);
  surfMean.assign(numFinal, 0.0);

  if(doMergeTree_) {
    this->computeMergeTree<dataType, triangulationType>(
      triangulation, finalTraj, surfMin, surfMax, surfMean, vertexTrajPerFrame);
  } else {
    vertexTrajPerFrame.clear();
  }

  this->printMsg("Post-processing complete", 1.0, timer.getElapsedTime(),
                 this->threadNumber_);
  return 1;
}
