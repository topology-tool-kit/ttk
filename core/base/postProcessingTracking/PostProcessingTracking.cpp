#include <PostProcessingTracking.h>

ttk::PostProcessingTracking::PostProcessingTracking() {
  this->setDebugMsgPrefix("PostProcessingTracking");
}

int ttk::PostProcessingTracking::correctTrajectory(
  const std::vector<std::vector<int>> &trajTime,
  const std::vector<std::vector<int>> &trajVertexId,
  const std::vector<std::vector<double>> &coordsX,
  const std::vector<std::vector<double>> &coordsY,
  const std::vector<int> &trajCriticalType,
  std::vector<LinearTrajectory> &linearTraj,
  std::vector<LinearTrajectory> &outputTraj,
  std::vector<FuseRecord> &fuseRecords) {

  ttk::Timer timer;
  const int numTraj = static_cast<int>(trajTime.size());
  const bool useTypeConstraint
    = (static_cast<int>(trajCriticalType.size()) == numTraj);
  this->printMsg("Linearization and chaining (" + std::to_string(numTraj)
                 + " input trajectories, linearize="
                 + std::to_string(doLinearize_) + ", fuse="
                 + std::to_string(doFusion_) + ", linearizeFuse="
                 + std::to_string(doLinearizeFuse_) + ", same-type="
                 + std::to_string(useTypeConstraint) + ")");

  auto dirDot = [&](int i, int j, const std::vector<double> &mDx,
                    const std::vector<double> &mDy,
                    const std::vector<double> &mDz) -> double {
    return mDx[i] * mDx[j] + mDy[i] * mDy[j] + mDz[i] * mDz[j];
  };

  auto temporalOk = [&](int sFrame, int eFrame) -> bool {
    return (sFrame - eFrame > minFrameDist_)
           && (sFrame - eFrame < maxFrameDist_);
  };

  auto dist2AtStartFrame = [&](const LinearTrajectory &coefI,
                               const LinearTrajectory &coefJ,
                               int sFrame) -> double {
    const double xTh = coefI.evalX(sFrame);
    const double yTh = coefI.evalY(sFrame);
    const double xJ = coefJ.evalX(sFrame);
    const double yJ = coefJ.evalY(sFrame);
    const double dx = xJ - xTh, dy = yJ - yTh;
    return dx * dx + dy * dy;
  };

  auto buildSamplesForChain = [&](const std::vector<FuseRecord> &chain,
                                  std::vector<int> &T,
                                  std::vector<double> &X,
                                  std::vector<double> &Y) {
    const int capacity = static_cast<int>(chain.size()) * 2 + 2;
    T.reserve(capacity);
    X.reserve(capacity);
    Y.reserve(capacity);
    for(const auto &r : chain) {
      const std::vector<int> T2{trajTime[r.i].front(), r.endFrame};
      const auto &cI = linearTraj[r.i];
      for(const int t : T2) {
        X.push_back(cI.evalX(t));
        Y.push_back(cI.evalY(t));
        T.push_back(t);
      }
    }
    const FuseRecord &r = chain.back();
    const auto &cJ = linearTraj[r.j];
    X.push_back(cJ.evalX(r.startFrame));
    X.push_back(cJ.evalX(trajTime[r.j].back()));
    Y.push_back(cJ.evalY(r.startFrame));
    Y.push_back(cJ.evalY(trajTime[r.j].back()));
    T.push_back(r.startFrame);
    T.push_back(trajTime[r.j].back());
  };

  auto fitLineCoefForChain
    = [&](const std::vector<FuseRecord> &chain) -> LinearTrajectory {
    std::vector<int> T;
    std::vector<double> X, Y;
    buildSamplesForChain(chain, T, X, Y);
    LinearTrajectory lineCoef;
#ifdef TTK_ENABLE_EIGEN
    linearRegression(T, X, Y, lineCoef);
#else
    lineCoef = linearTraj[chain.front().i];
#endif
    lineCoef.startFrame = trajTime[chain[0].i].front();
    lineCoef.endFrame = trajTime[chain.back().j].back();
    return lineCoef;
  };

  linearTraj.assign(numTraj, LinearTrajectory{});
  for(int i = 0; i < numTraj; ++i) {
    if(trajTime[i].empty())
      continue;
    linearTraj[i].startFrame = trajTime[i].front();
    linearTraj[i].endFrame = trajTime[i].back();
    if(doLinearize_) {
#ifdef TTK_ENABLE_EIGEN
      linearRegression(trajTime[i], coordsX[i], coordsY[i], linearTraj[i]);
#else
      this->printWrn("Eigen unavailable: skipping linear regression");
#endif
    } else {
      const int tA = trajTime[i].front();
      const int tB = trajTime[i].back();
      const double xA = coordsX[i].front();
      const double xB = coordsX[i].back();
      const double yA = coordsY[i].front();
      const double yB = coordsY[i].back();
      if(tB != tA) {
        linearTraj[i].ax = (xB - xA) / static_cast<double>(tB - tA);
        linearTraj[i].ay = (yB - yA) / static_cast<double>(tB - tA);
        linearTraj[i].bx = xA - linearTraj[i].ax * tA;
        linearTraj[i].by = yA - linearTraj[i].ay * tA;
      } else {
        linearTraj[i].ax = 0.0;
        linearTraj[i].bx = xA;
        linearTraj[i].ay = 0.0;
        linearTraj[i].by = yA;
      }
    }
  }

  outputTraj.clear();
  fuseRecords.clear();

  if(!doFusion_) {
    outputTraj.reserve(numTraj);
    for(int i = 0; i < numTraj; ++i) {
      if(trajTime[i].empty())
        continue;
      LinearTrajectory lt = linearTraj[i];
      lt.originalTrajId = i;
      lt.segmentKind = 0;
      lt.criticalPoints.reserve(trajTime[i].size());
      for(size_t k = 0; k < trajTime[i].size(); ++k) {
        lt.criticalPoints.emplace_back(
          trajTime[i][k],
          static_cast<ttk::SimplexId>(trajVertexId[i][k]));
      }
      linearTraj[i].finalChainId = static_cast<int>(outputTraj.size());
      lt.finalChainId = linearTraj[i].finalChainId;
      outputTraj.push_back(std::move(lt));
    }
    this->printMsg("Linearization only (" + std::to_string(outputTraj.size())
                     + " output chains)",
                   1.0, timer.getElapsedTime(), this->threadNumber_);
    return 1;
  }

  std::vector<double> meanDx(numTraj), meanDy(numTraj), meanDz(numTraj);
  computeMeanUnitDirectionLinear(linearTraj, meanDx, meanDy, meanDz);

  fuseRecords.reserve(numTraj);
  std::vector<char> usedAsStart(numTraj, false), usedAsEnd(numTraj, false);

  const double similarityThreshold = cosCol_;
  const double maxLinkDist2 = maxRadius_;

  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || trajTime[i].empty()) continue;

    const int endFrame = trajTime[i].back();

    int    bestJ     = -1;
    double bestScore = std::numeric_limits<double>::infinity();

    const double dist2Denom  = (maxLinkDist2 > 0.0) ? maxLinkDist2 : 1.0;
    const double dotDenom    = (1.0 - similarityThreshold > 1e-12)
                               ? (1.0 - similarityThreshold) : 1.0;
    const double timeDenom   = (maxFrameDist_ > 0) ? static_cast<double>(maxFrameDist_) : 1.0;

    for(int j = 0; j < numTraj; ++j) {
      if(usedAsEnd[j] || j == i || trajTime[j].empty()) continue;

      if(useTypeConstraint
         && trajCriticalType[i] != trajCriticalType[j]) continue;

      const int startFrame = trajTime[j].front();

      const double dist2 = dist2AtStartFrame(linearTraj[i], linearTraj[j], startFrame);
      if(dist2 > maxLinkDist2) continue;

      const double dot = dirDot(i, j, meanDx, meanDy, meanDz);
      if(dot < similarityThreshold) continue;

      if(!temporalOk(startFrame, endFrame)) continue;

      const double penDist2 = dist2 / dist2Denom;
      const double penDot   = (1.0 - dot) / dotDenom;
      const double penTime  = static_cast<double>(std::abs(endFrame - startFrame))
                              / timeDenom;

      const double score = penDist2 + penDot + penTime;

      if(score < bestScore) {
        bestScore = score;
        bestJ     = j;
      }
    }
    if(bestJ >= 0) {
      fuseRecords.push_back(
        {i, bestJ, trajTime[i].back(), trajTime[bestJ].front(), -1});
      usedAsStart[i] = true;
      usedAsEnd[bestJ] = true;
    }
  }

  outputTraj.reserve(numTraj);
  std::vector<bool> used(fuseRecords.size(), false);

  for(size_t idx1 = 0; idx1 < fuseRecords.size(); ++idx1) {
    if(used[idx1])
      continue;

    auto &r1 = fuseRecords[idx1];
    const int finalId = static_cast<int>(outputTraj.size());

    std::vector<FuseRecord> chain{r1};
    r1.finalContrib = finalId;
    linearTraj[r1.i].finalChainId = finalId;
    linearTraj[r1.j].finalChainId = finalId;
    used[idx1] = true;

    // Prepend 
    bool prepended = true;
    while(prepended) {
      prepended = false;
      for(size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
        if(used[idx2])
          continue;
        auto &r2 = fuseRecords[idx2];
        if(r2.j == chain.front().i) {
          chain.insert(chain.begin(), r2);
          r2.finalContrib = finalId;
          linearTraj[r2.i].finalChainId = finalId;
          linearTraj[r2.j].finalChainId = finalId;
          used[idx2] = true;
          prepended = true;
          break;
        }
      }
    }

    // Extend 
    bool extended = true;
    while(extended) {
      extended = false;
      for(size_t idx2 = 0; idx2 < fuseRecords.size(); ++idx2) {
        if(used[idx2])
          continue;
        auto &r2 = fuseRecords[idx2];
        if(chain.back().j == r2.i) {
          chain.push_back(r2);
          r2.finalContrib = finalId;
          linearTraj[r2.i].finalChainId = finalId;
          linearTraj[r2.j].finalChainId = finalId;
          used[idx2] = true;
          extended = true;
          break;
        }
      }
    }

    if(doLinearizeFuse_) {
      LinearTrajectory lineCoef = fitLineCoefForChain(chain);
      lineCoef.finalChainId = finalId;
      lineCoef.originalTrajId = -1;
      lineCoef.segmentKind = 2;

      const int firstTraj = chain[0].i;
      for(size_t k = 0; k < trajTime[firstTraj].size(); ++k) {
        lineCoef.criticalPoints.emplace_back(
          trajTime[firstTraj][k],
          static_cast<ttk::SimplexId>(trajVertexId[firstTraj][k]));
      }
      for(const auto &rec : chain) {
        const int tj = rec.j;
        for(size_t k = 0; k < trajTime[tj].size(); ++k) {
          lineCoef.criticalPoints.emplace_back(
            trajTime[tj][k],
            static_cast<ttk::SimplexId>(trajVertexId[tj][k]));
        }
      }

      outputTraj.push_back(std::move(lineCoef));
    } else {
      const int firstTraj = chain[0].i;
      {
        LinearTrajectory seg = linearTraj[firstTraj];
        seg.startFrame = trajTime[firstTraj].front();
        seg.endFrame = trajTime[firstTraj].back();
        seg.finalChainId = finalId;
        seg.originalTrajId = firstTraj;
        seg.segmentKind = 0;
        seg.criticalPoints.reserve(trajTime[firstTraj].size());
        for(size_t k = 0; k < trajTime[firstTraj].size(); ++k) {
          seg.criticalPoints.emplace_back(
            trajTime[firstTraj][k],
            static_cast<ttk::SimplexId>(trajVertexId[firstTraj][k]));
        }
        outputTraj.push_back(std::move(seg));
      }

      for(const auto &rec : chain) {
        const int iSeg = rec.i;
        const int jSeg = rec.j;
        const int tEnd = trajTime[iSeg].back();
        const int tStart = trajTime[jSeg].front();
        const auto &cI = linearTraj[iSeg];
        const auto &cJ = linearTraj[jSeg];

        LinearTrajectory junction;
        const double xA = cI.evalX(tEnd);
        const double yA = cI.evalY(tEnd);
        const double xB = cJ.evalX(tStart);
        const double yB = cJ.evalY(tStart);
        if(tStart != tEnd) {
          const double dt = static_cast<double>(tStart - tEnd);
          junction.ax = (xB - xA) / dt;
          junction.ay = (yB - yA) / dt;
          junction.bx = xA - junction.ax * static_cast<double>(tEnd);
          junction.by = yA - junction.ay * static_cast<double>(tEnd);
        } else {
          junction.ax = 0.0;
          junction.ay = 0.0;
          junction.bx = xA;
          junction.by = yA;
        }
        junction.startFrame = tEnd;
        junction.endFrame = tStart;
        junction.finalChainId = finalId;
        junction.originalTrajId = -1;
        junction.segmentKind = 1;
        outputTraj.push_back(std::move(junction));

        LinearTrajectory seg = linearTraj[jSeg];
        seg.startFrame = trajTime[jSeg].front();
        seg.endFrame = trajTime[jSeg].back();
        seg.finalChainId = finalId;
        seg.originalTrajId = jSeg;
        seg.segmentKind = 0;
        seg.criticalPoints.reserve(trajTime[jSeg].size());
        for(size_t k = 0; k < trajTime[jSeg].size(); ++k) {
          seg.criticalPoints.emplace_back(
            trajTime[jSeg][k],
            static_cast<ttk::SimplexId>(trajVertexId[jSeg][k]));
        }
        outputTraj.push_back(std::move(seg));
      }
    }
  }

  for(int i = 0; i < numTraj; ++i) {
    if(usedAsStart[i] || usedAsEnd[i] || trajTime[i].empty())
      continue;
    LinearTrajectory lineCoef = linearTraj[i];
    lineCoef.startFrame = trajTime[i].front();
    lineCoef.endFrame = trajTime[i].back();
    lineCoef.originalTrajId = i;
    lineCoef.segmentKind = 0;
    lineCoef.criticalPoints.reserve(trajTime[i].size());
    for(size_t k = 0; k < trajTime[i].size(); ++k) {
      lineCoef.criticalPoints.emplace_back(
        trajTime[i][k],
        static_cast<ttk::SimplexId>(trajVertexId[i][k]));
    }
    linearTraj[i].finalChainId = static_cast<int>(outputTraj.size());
    lineCoef.finalChainId = linearTraj[i].finalChainId;
    outputTraj.push_back(std::move(lineCoef));
  }

  for(auto &c : outputTraj) {
    if(c.endFrame < c.startFrame)
      std::swap(c.startFrame, c.endFrame);
  }

  this->printMsg("Linearization and chaining ("
                   + std::to_string(outputTraj.size()) + " output chains)",
                 1.0, timer.getElapsedTime(), this->threadNumber_);
  return 1;
}
