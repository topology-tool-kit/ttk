/// \ingroup base
/// \class ttk::MergeTreeTemporalReductionDecoding
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// This module defines the %MergeTreeTemporalReductionDecoding class that
/// computes the reconstruction of a reduced sequence of merge trees.
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
   * The MergeTreeTemporalReductionDecoding class provides methods to compute
   * the reconstruction of a reduced sequence of merge trees.
   */
  class MergeTreeTemporalReductionDecoding : virtual public Debug,
                                             public MergeTreeBase {
  protected:
    std::vector<double> finalDistances_, distancesToKeyFrames_;

  public:
    MergeTreeTemporalReductionDecoding();

    template <class dataType>
    dataType computeDistance(
      ftm::MergeTree<dataType> &mTree1,
      ftm::MergeTree<dataType> &mTree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching) {
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

      dataType distance
        = mergeTreeDistance.execute<dataType>(mTree1, mTree2, matching);

      return distance;
    }

    template <class dataType>
    dataType computeDistance(ftm::MergeTree<dataType> &mTree1,
                             ftm::MergeTree<dataType> &mTree2) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      return computeDistance<dataType>(mTree1, mTree2, matching);
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

    template <class dataType>
    void execute(
      std::vector<ftm::MergeTree<dataType>> &mTrees,
      std::vector<std::tuple<double, int, int, int, int>> &coefs,
      std::vector<ftm::MergeTree<dataType>> &allMT,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &allMatching) {
      Timer t_tempSub;

      // --- Preprocessing
      treesNodeCorr_ = std::vector<std::vector<int>>(mTrees.size());
      for(unsigned int i = 0; i < mTrees.size(); ++i) {
        preprocessingPipeline<dataType>(
          mTrees[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[i]);
      }
      printTreesStats<dataType>(mTrees);

      // --- Execute
      distancesToKeyFrames_ = std::vector<double>(coefs.size() * 2);
      int index = 0;
      size_t cpt = 0;
      while(cpt < coefs.size()) {
        while(cpt < coefs.size() and std::get<2>(coefs[cpt]) <= index) {
          double alpha = std::get<0>(coefs[cpt]);
          int index1 = std::get<1>(coefs[cpt]);
          int index2 = std::get<2>(coefs[cpt]);
          ftm::MergeTree<dataType> tree = computeBarycenter<dataType>(
            mTrees[index1], mTrees[index2], alpha);
          allMT.push_back(tree);
          distancesToKeyFrames_[cpt * 2]
            = computeDistance<dataType>(mTrees[index1], tree);
          distancesToKeyFrames_[cpt * 2 + 1]
            = computeDistance<dataType>(tree, mTrees[index2]);
          ++cpt;
        }
        allMT.push_back(mTrees[index]);
        ++index;
      }

      allMatching = std::vector<
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>(allMT.size()
                                                                   - 1);
      finalDistances_ = std::vector<double>(allMT.size() - 1);
      for(unsigned int i = 0; i < allMT.size() - 1; ++i)
        finalDistances_[i]
          = computeDistance<dataType>(allMT[i], allMT[i + 1], allMatching[i]);

      // --- Postprocessing
      for(unsigned int i = 0; i < allMT.size(); ++i)
        postprocessingPipeline<dataType>(&(allMT[i].tree));
      for(unsigned int i = 0; i < mTrees.size(); ++i)
        postprocessingPipeline<dataType>(&(mTrees[i].tree));

      // --- Print results
      std::stringstream ss, ss2, ss3;
      ss << "input size    = " << mTrees.size();
      printMsg(ss.str());
      ss2 << "output size   = " << allMT.size();
      printMsg(ss2.str());
      ss3 << "reconstructed : " << allMT.size() - mTrees.size();
      printMsg(ss3.str());
      printMsg("Decoding", 1, t_tempSub.getElapsedTime(), this->threadNumber_);
    }

  }; // MergeTreeTemporalReductionDecoding class

} // namespace ttk
