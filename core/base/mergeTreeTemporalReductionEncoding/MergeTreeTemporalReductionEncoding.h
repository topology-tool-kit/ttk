/// \ingroup base
/// \class ttk::MergeTreeTemporalReductionEncoding
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// This module defines the %MergeTreeTemporalReductionEncoding class that
/// computes a temporal reduction of a sequence of merge trees.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021

#pragma once

// ttk common includes
#include <Debug.h>

#include <FTMTreeUtils.h>
#include <MergeTreeBarycenter.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>

namespace ttk {

  /**
   * The MergeTreeTemporalReductionEncoding class provides methods to compute
   * a temporal reduction of a sequence of merge trees.
   */
  class MergeTreeTemporalReductionEncoding : virtual public Debug,
                                             public MergeTreeBase {
  protected:
    double removalPercentage_ = 50.;
    bool useL2Distance_ = false;
    std::vector<std::vector<double>> fieldL2_;
    bool useCustomTimeVariable_ = false;
    std::vector<double> timeVariable_;

  public:
    MergeTreeTemporalReductionEncoding();

    void setRemovalPercentage(double rs) {
      removalPercentage_ = rs;
    }

    void setUseL2Distance(bool useL2) {
      useL2Distance_ = useL2;
    }

    template <class dataType>
    dataType computeL2Distance(std::vector<dataType> &img1,
                               std::vector<dataType> &img2,
                               bool emptyFieldDistance = false) {
      size_t noPoints = img1.size();

      std::vector<dataType> secondField = img2;
      if(emptyFieldDistance)
        secondField = std::vector<dataType>(noPoints, 0);

      dataType distance = 0;

      for(size_t i = 0; i < noPoints; ++i)
        distance += std::pow((img1[i] - secondField[i]), 2);

      distance = std::sqrt(distance);

      return distance;
    }

    template <class dataType>
    std::vector<dataType> computeL2Barycenter(std::vector<dataType> &img1,
                                              std::vector<dataType> &img2,
                                              double alpha) {

      size_t noPoints = img1.size();

      std::vector<dataType> barycenter(noPoints);
      for(size_t i = 0; i < noPoints; ++i)
        barycenter[i] = alpha * img1[i] * (1 - alpha) * img2[i];

      return barycenter;
    }

    template <class dataType>
    dataType computeDistance(ftm::MergeTree<dataType> &mTree1,
                             ftm::MergeTree<dataType> &mTree2,
                             bool emptyTreeDistance = false) {
      MergeTreeDistance mergeTreeDistance;
      mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
      mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
      mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
      mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
      mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
      mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
      mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
      mergeTreeDistance.setProgressiveComputation(progressiveComputation_);
      mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
      mergeTreeDistance.setParallelize(parallelize_);
      mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
      mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeDistance.setNormalizedWassersteinReg(normalizedWassersteinReg_);
      mergeTreeDistance.setRescaledWasserstein(rescaledWasserstein_);
      mergeTreeDistance.setKeepSubtree(keepSubtree_);
      mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
      mergeTreeDistance.setThreadNumber(this->threadNumber_);
      mergeTreeDistance.setDistanceSquared(true); // squared root
      mergeTreeDistance.setDebugLevel(2);
      mergeTreeDistance.setPreprocess(false);
      mergeTreeDistance.setPostprocess(false);
      // mergeTreeDistance.setIsCalled(true);
      mergeTreeDistance.setOnlyEmptyTreeDistance(emptyTreeDistance);

      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      dataType distance
        = mergeTreeDistance.execute<dataType>(mTree1, mTree2, matching);

      return distance;
    }

    template <class dataType>
    ftm::MergeTree<dataType> computeBarycenter(ftm::MergeTree<dataType> &mTree1,
                                               ftm::MergeTree<dataType> &mTree2,
                                               double alpha) {
      MergeTreeBarycenter mergeTreeBarycenter;
      mergeTreeBarycenter.setAssignmentSolver(assignmentSolverID_);
      mergeTreeBarycenter.setEpsilonTree1(epsilonTree1_);
      mergeTreeBarycenter.setEpsilonTree2(epsilonTree2_);
      mergeTreeBarycenter.setEpsilon2Tree1(epsilon2Tree1_);
      mergeTreeBarycenter.setEpsilon2Tree2(epsilon2Tree2_);
      mergeTreeBarycenter.setEpsilon3Tree1(epsilon3Tree1_);
      mergeTreeBarycenter.setEpsilon3Tree2(epsilon3Tree2_);
      mergeTreeBarycenter.setProgressiveComputation(progressiveComputation_);
      mergeTreeBarycenter.setBranchDecomposition(branchDecomposition_);
      mergeTreeBarycenter.setParallelize(parallelize_);
      mergeTreeBarycenter.setPersistenceThreshold(persistenceThreshold_);
      mergeTreeBarycenter.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeBarycenter.setNormalizedWassersteinReg(
        normalizedWassersteinReg_);
      mergeTreeBarycenter.setRescaledWasserstein(rescaledWasserstein_);
      mergeTreeBarycenter.setKeepSubtree(keepSubtree_);
      mergeTreeBarycenter.setUseMinMaxPair(useMinMaxPair_);
      mergeTreeBarycenter.setThreadNumber(this->threadNumber_);
      mergeTreeBarycenter.setAlpha(alpha);
      mergeTreeBarycenter.setDebugLevel(2);
      mergeTreeBarycenter.setPreprocess(false);
      mergeTreeBarycenter.setPostprocess(false);
      // mergeTreeBarycenter.setIsCalled(true);

      std::vector<ftm::MergeTree<dataType>> intermediateTrees;
      intermediateTrees.push_back(mTree1);
      intermediateTrees.push_back(mTree2);
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        outputMatchingBarycenter(2);
      ftm::MergeTree<dataType> barycenter;
      mergeTreeBarycenter.execute<dataType>(
        intermediateTrees, outputMatchingBarycenter, barycenter);
      return barycenter;
    }

    double computeAlpha(int index1, int middleIndex, int index2) {
      index1 = timeVariable_[index1];
      middleIndex = timeVariable_[middleIndex];
      index2 = timeVariable_[index2];
      return 1 - ((double)middleIndex - index1) / (index2 - index1);
    }

    template <class dataType>
    void
      temporalSubsampling(std::vector<ftm::MergeTree<dataType>> &mTrees,
                          std::vector<int> &removed,
                          std::vector<ftm::MergeTree<dataType>> &barycenters,
                          std::vector<std::vector<dataType>> &barycentersL2) {
      std::vector<bool> treeRemoved(mTrees.size(), false);

      int toRemoved = mTrees.size() * removalPercentage_ / 100.;
      toRemoved = std::min(toRemoved, (int)(mTrees.size() - 3));

      std::vector<std::vector<dataType>> images(fieldL2_.size());
      for(size_t i = 0; i < fieldL2_.size(); ++i)
        for(size_t j = 0; j < fieldL2_[i].size(); ++j)
          images[i].push_back(static_cast<dataType>(fieldL2_[i][j]));

      for(int iter = 0; iter < toRemoved; ++iter) {
        dataType bestCost = std::numeric_limits<dataType>::max();
        int bestMiddleIndex = -1;
        ftm::MergeTree<dataType> bestBarycenter;
        std::vector<std::tuple<ftm::MergeTree<dataType>, int>>
          bestBarycentersOnPath;
        std::vector<dataType> bestBarycenterL2;
        std::vector<std::tuple<std::vector<dataType>, int>>
          bestBarycentersL2OnPath;

        // Compute barycenter for each pair of trees
        printMsg("Compute barycenter for each pair of trees",
                 debug::Priority::VERBOSE);
        unsigned int index1 = 0, index2 = 0;
        while(index2 != mTrees.size() - 1) {

          // Get index in the middle
          int middleIndex = index1 + 1;
          while(treeRemoved[middleIndex])
            ++middleIndex;

          // Get second index
          index2 = middleIndex + 1;
          while(treeRemoved[index2])
            ++index2;

          // Compute barycenter
          printMsg("Compute barycenter", debug::Priority::VERBOSE);
          double alpha = computeAlpha(index1, middleIndex, index2);
          ftm::MergeTree<dataType> barycenter;
          std::vector<dataType> barycenterL2;
          if(not useL2Distance_)
            barycenter = computeBarycenter<dataType>(
              mTrees[index1], mTrees[index2], alpha);
          else
            barycenterL2 = computeL2Barycenter<dataType>(
              images[index1], images[index2], alpha);

          // - Compute cost
          // Compute distance with middleIndex
          printMsg(
            "Compute distance with middleIndex", debug::Priority::VERBOSE);
          dataType cost;
          if(not useL2Distance_)
            cost = computeDistance<dataType>(barycenter, mTrees[middleIndex]);
          else
            cost
              = computeL2Distance<dataType>(barycenterL2, images[middleIndex]);

          // Compute distances of previously removed trees on the path
          printMsg("Compute distances of previously removed trees",
                   debug::Priority::VERBOSE);
          std::vector<std::tuple<ftm::MergeTree<dataType>, int>>
            barycentersOnPath;
          std::vector<std::tuple<std::vector<dataType>, int>>
            barycentersL2OnPath;
          for(unsigned int i = 0; i < 2; ++i) {
            int toReach = (i == 0 ? index1 : index2);
            int offset = (i == 0 ? -1 : 1);
            int tIndex = middleIndex + offset;
            while(tIndex != toReach) {

              // Compute barycenter
              double alphaT = computeAlpha(index1, tIndex, index2);
              ftm::MergeTree<dataType> barycenterP;
              std::vector<dataType> barycenterPL2;
              if(not useL2Distance_)
                barycenterP = computeBarycenter<dataType>(
                  mTrees[index1], mTrees[index2], alphaT);
              else
                barycenterPL2 = computeL2Barycenter<dataType>(
                  images[index1], images[index2], alphaT);

              // Compute distance
              dataType costP;
              if(not useL2Distance_)
                costP = computeDistance<dataType>(barycenterP, mTrees[tIndex]);
              else
                costP
                  = computeL2Distance<dataType>(barycenterPL2, images[tIndex]);

              // Save results
              if(not useL2Distance_)
                barycentersOnPath.push_back(
                  std::make_tuple(barycenterP, tIndex));
              else
                barycentersL2OnPath.push_back(
                  std::make_tuple(barycenterPL2, tIndex));
              cost += costP;
              tIndex += offset;
            }
          }

          if(cost < bestCost) {
            bestCost = cost;
            bestMiddleIndex = middleIndex;
            if(not useL2Distance_) {
              bestBarycenter = barycenter;
              bestBarycentersOnPath = barycentersOnPath;
            } else {
              bestBarycenterL2 = barycenterL2;
              bestBarycentersL2OnPath = barycentersL2OnPath;
            }
          }

          // Go to the next index
          index1 = middleIndex;
        }

        // Removed the tree with the lowest cost
        printMsg(
          "Removed the tree with the lowest cost", debug::Priority::VERBOSE);
        removed.push_back(bestMiddleIndex);
        treeRemoved[bestMiddleIndex] = true;
        if(not useL2Distance_) {
          barycenters[bestMiddleIndex] = bestBarycenter;
          for(auto &tup : bestBarycentersOnPath)
            barycenters[std::get<1>(tup)] = std::get<0>(tup);
        } else {
          barycentersL2[bestMiddleIndex] = bestBarycenterL2;
          for(auto &tup : bestBarycentersL2OnPath)
            barycentersL2[std::get<1>(tup)] = std::get<0>(tup);
        }
      }
    }

    template <class dataType>
    std::vector<int> execute(std::vector<ftm::MergeTree<dataType>> &mTrees,
                             std::vector<double> &emptyTreeDistances,
                             std::vector<ftm::MergeTree<dataType>> &allMT) {
      Timer t_tempSub;

      // --- Preprocessing
      if(not useL2Distance_) {
        treesNodeCorr_ = std::vector<std::vector<int>>(mTrees.size());
        for(unsigned int i = 0; i < mTrees.size(); ++i) {
          preprocessingPipeline<dataType>(mTrees[i], epsilonTree2_,
                                          epsilon2Tree2_, epsilon3Tree2_,
                                          branchDecomposition_, useMinMaxPair_,
                                          cleanTree_, treesNodeCorr_[i]);
        }
        printTreesStats<dataType>(mTrees);
      }

      // --- Execute
      std::vector<ftm::MergeTree<dataType>> barycenters(mTrees.size());
      std::vector<std::vector<dataType>> barycentersL2(mTrees.size());
      std::vector<int> removed;
      if(not useCustomTimeVariable_) {
        timeVariable_.clear();
        for(size_t i = 0; i < mTrees.size(); ++i)
          timeVariable_.push_back(i);
      }
      temporalSubsampling<dataType>(
        mTrees, removed, barycenters, barycentersL2);

      // --- Concatenate all trees/L2Images
      std::vector<std::vector<dataType>> images(fieldL2_.size());
      for(size_t i = 0; i < fieldL2_.size(); ++i)
        for(size_t j = 0; j < fieldL2_[i].size(); ++j)
          images[i].push_back(static_cast<dataType>(fieldL2_[i][j]));

      for(auto &mt : mTrees)
        allMT.push_back(mt);
      std::vector<bool> removedB(mTrees.size(), false);
      for(auto r : removed)
        removedB[r] = true;
      for(unsigned int i = 0; i < barycenters.size(); ++i)
        if(removedB[i]) {
          if(not useL2Distance_)
            allMT.push_back(barycenters[i]);
          else
            images.push_back(barycentersL2[i]);
        }

      // --- Compute empty tree distances
      unsigned int distMatSize
        = (not useL2Distance_ ? allMT.size() : images.size());
      for(unsigned int i = 0; i < distMatSize; ++i) {
        dataType distance;
        if(not useL2Distance_)
          distance = computeDistance<dataType>(allMT[i], allMT[i], true);
        else
          distance = computeL2Distance<dataType>(images[i], images[i], true);
        emptyTreeDistances.push_back(distance);
      }

      // --- Postprocessing
      if(not useL2Distance_) {
        for(unsigned int i = 0; i < allMT.size(); ++i)
          postprocessingPipeline<dataType>(&(allMT[i].tree));
        for(unsigned int i = 0; i < mTrees.size(); ++i)
          postprocessingPipeline<dataType>(&(mTrees[i].tree));
      }

      // --- Print results
      std::stringstream ss, ss2, ss3;
      ss << "input size    = " << mTrees.size();
      printMsg(ss.str());
      ss2 << "output size   = "
          << mTrees.size() - (distMatSize - mTrees.size());
      printMsg(ss2.str());
      ss3 << "removed       : ";
      for(unsigned int i = 0; i < removed.size(); ++i) {
        auto r = removed[i];
        ss3 << r;
        if(i < removed.size() - 1)
          ss3 << ", ";
      }
      printMsg(ss3.str());

      sort(removed.begin(), removed.end());

      printMsg("Encoding", 1, t_tempSub.getElapsedTime(), this->threadNumber_);

      return removed;
    }

  }; // MergeTreeTemporalReductionEncoding class

} // namespace ttk
