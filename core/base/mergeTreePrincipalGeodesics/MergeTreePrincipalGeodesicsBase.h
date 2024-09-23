/// \ingroup base
/// \class ttk::MergeTreePrincipalGeodesicsBase
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022
///
/// \b Related \b publication: \n
/// "Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Jules Vidal, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2022

#pragma once

#include <MergeTreeAxesAlgorithmBase.h>
#include <cmath>

#define POINTS_BELOW_DIAG_TOLERANCE 1e-4

namespace ttk {

  class MergeTreePrincipalGeodesicsBase : virtual public Debug,
                                          public MergeTreeAxesAlgorithmBase {
  protected:
    // Filled by the algorithm
    std::vector<std::vector<std::vector<double>>> vS_, v2s_, trees2Vs_,
      trees2V2s_;
    std::vector<std::vector<double>> allTs_, allScaledTs_, allTreesTs_;
    std::vector<std::vector<double>> branchesCorrelationMatrix_,
      persCorrelationMatrix_;
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      baryMatchings_, baryMatchings2_;

  public:
    MergeTreePrincipalGeodesicsBase() {
      this->setDebugMsgPrefix(
        "MergeTreePrincipalGeodesicsBase"); // inherited from Debug: prefix will
                                            // be printed at the beginning of
                                            // every msg
    }

    //----------------------------------------------------------------------------
    // Matching / Distance
    //----------------------------------------------------------------------------
    template <class dataType>
    dataType computeReconstructionError(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &inputTrees,
      std::vector<std::vector<double *>> &vS,
      std::vector<std::vector<double *>> &v2s,
      size_t vSize,
      std::vector<std::vector<double>> &allTreesTs,
      std::vector<double> &reconstructionErrors,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      bool transposeVector = true) {
      reconstructionErrors.resize(inputTrees.size());
      matchings.resize(inputTrees.size());
      dataType reconstructionError = 0.0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_) reduction(+:reconstructionError)
#endif
      for(unsigned int i = 0; i < inputTrees.size(); ++i) {
        ftm::MergeTree<dataType> reconstructedTree;
        getMultiInterpolation<dataType>(barycenter, vS, v2s, vSize,
                                        allTreesTs[i], reconstructedTree,
                                        transposeVector);
        dataType error;
        computeOneDistance<dataType>(
          reconstructedTree, inputTrees[i], matchings[i], error, true);
        reconstructionError += error / inputTrees.size();
        reconstructionErrors[i] = error;
      }
      return reconstructionError;
    }

    template <class dataType>
    dataType computeReconstructionError(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &inputTrees,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<double>> &allTreesTs,
      std::vector<double> &reconstructionErrors) {
      std::vector<std::vector<double *>> pVS, pV2s;
      vectorOfVectorsToPointers(vS, pVS);
      vectorOfVectorsToPointers(v2s, pV2s);
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        matchings;
      return computeReconstructionError(barycenter, inputTrees, pVS, pV2s,
                                        vS[0].size(), allTreesTs,
                                        reconstructionErrors, matchings, false);
    }

    template <class dataType>
    dataType computeReconstructionError(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &inputTrees,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<double>> &allTreesTs) {
      std::vector<double> reconstructionErrors;
      return computeReconstructionError(
        barycenter, inputTrees, vS, v2s, allTreesTs, reconstructionErrors);
    }

    //----------------------------------------------------------------------------
    // Interpolation
    //----------------------------------------------------------------------------
    double getGeodesicVectorMiddle(std::vector<double> &v,
                                   std::vector<double> &v2) {
      std::vector<double> vProj, v2Proj;
      ttk::Geometry::addVectorsProjection(v, v2, vProj, v2Proj);

      double alpha = 0.0;
      int cptDivide = 0;
      for(unsigned int i = 0; i < vProj.size(); ++i) {
        if(std::abs(v2Proj[i]) < Geometry::pow(10.0, -DBL_DIG))
          continue;
        alpha += (vProj[i] / v2Proj[i]);
        ++cptDivide;
      }
      alpha /= cptDivide;
      double const m = alpha / (1 + alpha);
      return m;
    }

    double getAdjustedTMax(double tMin, double m) {
      if(std::isnan(m))
        return tMin;
      return (m - tMin) * (1 - m) / m + m;
    }

    double getAdjustedTMin(double tMax, double m) {
      if(std::isnan(m))
        return tMax;
      return (m - tMax) * m / (1 - m) + m;
    }

    void updateT(
      double newT, double m, double &tMin, double &tMax, bool updateTMin) {
      if(updateTMin) {
        tMin = std::max(newT, tMin);
        auto adjustedTMax = getAdjustedTMax(newT, m);
        if(adjustedTMax > tMin)
          tMax = std::min(adjustedTMax, tMax);
      } else {
        tMax = std::min(newT, tMax);
        auto adjustedTMin = getAdjustedTMin(newT, m);
        if(adjustedTMin < tMax)
          tMin = std::max(adjustedTMin, tMin);
      }
    }

    template <class dataType>
    bool adjustDiagonalT(dataType baryBirth,
                         dataType baryDeath,
                         ftm::idNode node,
                         std::vector<std::vector<double>> &vNew,
                         std::vector<std::vector<double>> &v2New,
                         double &tMin,
                         double &tMax) {
      bool shortener = false;

      // Compute V1 extremity
      auto newBirthV1 = baryBirth - vNew[node][0];
      auto newDeathV1 = baryDeath - vNew[node][1];

      // Compute V2 extremity
      auto newBirthV2 = baryBirth + v2New[node][0];
      auto newDeathV2 = baryDeath + v2New[node][1];

      // Adjust t and v
      if(newBirthV1 > newDeathV1 and newBirthV2 > newDeathV2) {
        shortener = true;

        double const divisor
          = (vNew[node][1] - vNew[node][0]) / (baryDeath - baryBirth);
        vNew[node][0] /= divisor;
        vNew[node][1] /= divisor;
        double const divisor2
          = (v2New[node][0] - v2New[node][1]) / (baryDeath - baryBirth);
        v2New[node][0] /= divisor2;
        v2New[node][1] /= divisor2;
      } else if(newBirthV1 > newDeathV1 or newBirthV2 > newDeathV2) {
        double const newT
          = (baryDeath - vNew[node][1] - baryBirth + vNew[node][0])
            / (vNew[node][0] + v2New[node][0]
               - (vNew[node][1] + v2New[node][1]));
        double const m = getGeodesicVectorMiddle(vNew[node], v2New[node]);
        if(newBirthV1 > newDeathV1)
          updateT(newT, m, tMin, tMax, true);
        if(newBirthV2 > newDeathV2)
          updateT(newT, m, tMin, tMax, false);
      }

      return shortener;
    }

    template <class dataType>
    bool adjustNestingT(ftm::MergeTree<dataType> &ttkNotUsed(barycenter),
                        dataType baryBirth,
                        dataType baryDeath,
                        ftm::idNode node,
                        std::vector<std::vector<double>> &vNew,
                        std::vector<std::vector<double>> &v2New,
                        double &tMin,
                        double &tMax) {
      bool shortener = false;

      dataType birthParent = 0.0;
      dataType deathParent = 1.0;

      auto interpolant = [&](int i, double t) {
        return -vNew[node][i] + t * (vNew[node][i] + v2New[node][i]);
      };

      // Compute V1 extremity
      auto newBirthV1 = baryBirth + interpolant(0, tMin);
      auto newDeathV1 = baryDeath + interpolant(1, tMin);
      auto parentBirthV1 = 0.0;
      auto parentDeathV1 = 1.0;

      // Compute V2 extremity
      auto newBirthV2 = baryBirth + interpolant(0, tMax);
      auto newDeathV2 = baryDeath + interpolant(1, tMax);
      auto parentBirthV2 = 0.0;
      auto parentDeathV2 = 1.0;

      //
      bool const extremOutPathIn
        = ((newDeathV1 > parentDeathV1 and not(newBirthV1 < parentBirthV1)
            and not(newDeathV2 > parentDeathV2) and newBirthV2 < parentBirthV2)
           or (not(newDeathV1 > parentDeathV1) and newBirthV1 < parentBirthV1
               and newDeathV2 > parentDeathV2
               and not(newBirthV2 < parentBirthV2)));

      // --- Adjust V
      // Both extremities outside and path outside
      if(((newDeathV1 > parentDeathV1 and newDeathV2 > parentDeathV2)
          or (newBirthV1 < parentBirthV1 and newBirthV2 < parentBirthV2))
         and not(extremOutPathIn)) {
        shortener = true;

        // Shorten V1
        double coef1 = 1.0;
        if(newDeathV1 > parentDeathV1 and newBirthV1 < parentBirthV1)
          if(newBirthV1 < newDeathV1)
            coef1 = (parentBirthV1 - baryBirth) / interpolant(0, tMin);
          else
            coef1 = (parentDeathV1 - baryDeath) / interpolant(1, tMin);
        else if(newBirthV1 < parentBirthV1)
          coef1 = (parentBirthV1 - baryBirth) / interpolant(0, tMin);
        else if(newDeathV1 > parentDeathV1)
          coef1 = (parentDeathV1 - baryDeath) / interpolant(1, tMin);
        vNew[node][0] *= coef1;
        vNew[node][1] *= coef1;

        // Shorten V2
        double coef2 = 1.0;
        if(newDeathV2 > parentDeathV2 and newBirthV2 < parentBirthV2)
          if(newBirthV2 < newDeathV2)
            coef2 = (parentBirthV2 - baryBirth) / interpolant(0, tMax);
          else
            coef2 = (parentDeathV2 - baryDeath) / interpolant(1, tMax);
        else if(newBirthV2 < parentBirthV2)
          coef2 = (parentBirthV2 - baryBirth) / interpolant(0, tMax);
        else if(newDeathV2 > parentDeathV2)
          coef2 = (parentDeathV2 - baryDeath) / interpolant(1, tMax);
        v2New[node][0] *= coef2;
        v2New[node][1] *= coef2;
      }

      // --- Adjust T
      // One extremity outside or both extremities outside but path inside
      if(((newDeathV1 > parentDeathV1 or newBirthV1 < parentBirthV1)
          != (newDeathV2 > parentDeathV2 or newBirthV2 < parentBirthV2))
         or extremOutPathIn) {
        double const m = getGeodesicVectorMiddle(vNew[node], v2New[node]);
        if(newDeathV1 > parentDeathV1 or newDeathV2 > parentDeathV2) {
          double const newT = (deathParent - baryDeath + vNew[node][1])
                              / (vNew[node][1] + v2New[node][1]);
          if(newDeathV1 > parentDeathV1)
            updateT(newT, m, tMin, tMax, true);
          if(newDeathV2 > parentDeathV2)
            updateT(newT, m, tMin, tMax, false);
        }
        if(newBirthV1 < parentBirthV1 or newBirthV2 < parentBirthV2) {
          double const newT = (birthParent - baryBirth + vNew[node][0])
                              / (vNew[node][0] + v2New[node][0]);
          if(newBirthV1 < parentBirthV1)
            updateT(newT, m, tMin, tMax, true);
          if(newBirthV2 < parentBirthV2)
            updateT(newT, m, tMin, tMax, false);
        }
      }

      return shortener;
    }

    template <class dataType>
    void getInterpolationVectorDebugMsg(dataType birth,
                                        dataType death,
                                        std::vector<std::vector<double>> &vNew,
                                        std::vector<std::vector<double>> &v2New,
                                        int i,
                                        double t,
                                        double tMin,
                                        double tMax,
                                        const std::string &msg,
                                        std::stringstream &ssT) {
      double const tNew = (t * (tMax - tMin) + tMin);
      std::streamsize const sSize = std::cout.precision();
      ssT << std::setprecision(12) << std::endl << msg << std::endl;
      ssT << "interBirth : "
          << birth - vNew[i][0] + tNew * (vNew[i][0] + v2New[i][0]) << " _ "
          << "interDeath : "
          << death - vNew[i][1] + tNew * (vNew[i][1] + v2New[i][1])
          << std::endl;
      ssT << "v : " << vNew[i][0] << " _ " << vNew[i][1] << std::endl;
      ssT << "v2 : " << v2New[i][0] << " _ " << v2New[i][1] << std::endl;
      ssT << "ts : " << tMin << " _ " << tMax << std::endl;
      ssT << "t : " << t << " _ tNew : " << tNew << std::endl;
      std::cout.precision(sSize);
    }

    template <class dataType>
    double getTNew(ftm::MergeTree<dataType> &barycenter,
                   std::vector<std::vector<double>> &v,
                   std::vector<std::vector<double>> &v2,
                   int i,
                   double t) {
      ftm::FTMTree_MT *baryTree = &(barycenter.tree);
      double tMin = 0.0, tMax = 1.0;

      // - Get node information
      auto birthDeath = baryTree->getBirthDeath<dataType>(i);
      auto birth = std::get<0>(birthDeath);
      auto death = std::get<1>(birthDeath);

      // - Get parent information
      auto parent = baryTree->getParentSafe(i);
      auto parentBirthDeath = baryTree->getBirthDeath<dataType>(parent);
      auto parentBirth = std::get<0>(parentBirthDeath);
      auto parentDeath = std::get<1>(parentBirthDeath);

      if(normalizedWasserstein_ and !baryTree->notNeedToNormalize(i)) {
        birth = (birth - parentBirth) / (parentDeath - parentBirth);
        death = (death - parentBirth) / (parentDeath - parentBirth);
      }

      // - Adjust Diagonal T
      bool diagonalShortener
        = adjustDiagonalT(birth, death, i, v, v2, tMin, tMax);

      // - Adjust Nesting T
      bool nestingShortener = false;
      if(normalizedWasserstein_ and !baryTree->notNeedToNormalize(i)) {
        nestingShortener
          = adjustNestingT(barycenter, birth, death, i, v, v2, tMin, tMax);
      }

      // - Compute new t
      double const tNew = t * (tMax - tMin) + tMin;

      if(diagonalShortener or nestingShortener)
        printWrn("[getTNew] shortener");

      return tNew;
    }

    template <class dataType>
    void getInterpolationVector(ftm::MergeTree<dataType> &barycenter,
                                std::vector<double *> &v,
                                std::vector<double *> &v2,
                                size_t vSize,
                                double t,
                                std::vector<dataType> &interpolationVectorT,
                                bool transposeVector) {
      // --- Init
      ftm::FTMTree_MT *baryTree = &(barycenter.tree);
      std::vector<dataType> scalarsVector;
      ttk::ftm::getTreeScalars<dataType>(barycenter, scalarsVector);
      std::vector<double> interpolationVector(scalarsVector.size());
      std::vector<std::vector<double>> vNew, v2New;
      pointersToVectors(v, vSize, vNew);
      pointersToVectors(v2, vSize, v2New);
      if(transposeVector) {
        std::vector<std::vector<double>> out;
        ttk::Geometry::transposeMatrix(vNew, out);
        vNew = out;
        ttk::Geometry::transposeMatrix(v2New, out);
        v2New = out;
      }

      // --- Tree traversal
      int cptBelowDiagonal = 0, cptNotNesting = 0;
      std::queue<ftm::idNode> queue;
      queue.emplace(baryTree->getRoot());
      int noNestingShortener = 0, noDiagonalShortener = 0;
      while(!queue.empty()) {
        ftm::idNode const i = queue.front();
        queue.pop();

        // - Get node information
        auto nodeBirthDeath = baryTree->getBirthDeathNode<dataType>(i);
        auto nodeBirth = std::get<0>(nodeBirthDeath);
        auto nodeDeath = std::get<1>(nodeBirthDeath);
        auto birth = scalarsVector[nodeBirth];
        auto death = scalarsVector[nodeDeath];

        // - Get parent information
        auto parent = baryTree->getParentSafe(i);
        auto parentNodeBirthDeath
          = baryTree->getBirthDeathNode<dataType>(parent);
        auto parentNodeBirth = std::get<0>(parentNodeBirthDeath);
        auto parentNodeDeath = std::get<1>(parentNodeBirthDeath);
        auto parentBirth = scalarsVector[parentNodeBirth];
        auto parentDeath = scalarsVector[parentNodeDeath];
        auto interParentBirth = interpolationVector[parentNodeBirth];
        auto interParentDeath = interpolationVector[parentNodeDeath];

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        std::stringstream ssT;
        std::streamsize const sSize2 = std::cout.precision();
        ssT << std::setprecision(12);
        ssT << "info" << std::endl;
        ssT << "birth death : " << birth << " ___ " << death << std::endl;
        ssT << "par. birth death : " << parentBirth << " ___ " << parentDeath
            << std::endl;
        ssT << "par. inter birth death : " << interParentBirth << " ___ "
            << interParentDeath << std::endl;

        if(normalizedWasserstein_ and !baryTree->notNeedToNormalize(i)) {
          birth = (birth - parentBirth) / (parentDeath - parentBirth);
          death = (death - parentBirth) / (parentDeath - parentBirth);
          ssT << "norm. birth death : " << birth << " ___ " << death
              << std::endl;
        }
        std::cout.precision(sSize2);
        // --------------------------------------------------------------------

        // - Adjust tMin and tMax
        double tMin = 0.0, tMax = 1.0;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        getInterpolationVectorDebugMsg<dataType>(birth, death, vNew, v2New, i,
                                                 t, tMin, tMax,
                                                 "before adjustDiagonalT", ssT);
        // --------------------------------------------------------------------

        // - Adjust Diagonal T
        bool const diagonalShortener
          = adjustDiagonalT(birth, death, i, vNew, v2New, tMin, tMax);
        noDiagonalShortener += diagonalShortener;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        getInterpolationVectorDebugMsg<dataType>(
          birth, death, vNew, v2New, i, t, tMin, tMax,
          "after adjustDiagonalT/before adjustNestingT", ssT);
        // --------------------------------------------------------------------

        // - Adjust Nesting T
        if(normalizedWasserstein_ and !baryTree->notNeedToNormalize(i)) {
          bool const nestingShortener = adjustNestingT(
            barycenter, birth, death, i, vNew, v2New, tMin, tMax);
          noNestingShortener += nestingShortener;

          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          getInterpolationVectorDebugMsg<dataType>(birth, death, vNew, v2New, i,
                                                   t, tMin, tMax,
                                                   "after adjustNestingT", ssT);
          // ------------------------------------------------------------------
        }

        // - Compute interpolation
        double const tNew = t * (tMax - tMin) + tMin;
        double interBirth
          = birth - vNew[i][0] + tNew * (vNew[i][0] + v2New[i][0]);
        double interDeath
          = death - vNew[i][1] + tNew * (vNew[i][1] + v2New[i][1]);
        if(normalizedWasserstein_ and !baryTree->notNeedToNormalize(i)) {
          double const coef = (interParentDeath - interParentBirth);
          interBirth *= coef;
          interBirth += interParentBirth;
          interDeath *= coef;
          interDeath += interParentBirth;
        }
        interpolationVector[nodeBirth] = interBirth;
        interpolationVector[nodeDeath] = interDeath;

        // - Verify NaN
        if(std::isnan(interpolationVector[nodeBirth])
           or std::isnan(interpolationVector[nodeDeath])) {
          printMsg(ssT.str());
          printErr("NOT A NUMBER");
        }
        // - Verify points below diagonal
        double const eps = POINTS_BELOW_DIAG_TOLERANCE;
        if(interpolationVector[nodeBirth] > interpolationVector[nodeDeath]
           and interpolationVector[nodeBirth] - interpolationVector[nodeDeath]
                 < eps) {
          interpolationVector[nodeBirth] = interpolationVector[nodeDeath];
        }
        if(interpolationVector[nodeBirth] > interpolationVector[nodeDeath]) {
          ++cptBelowDiagonal;
          printMsg(ssT.str());
          std::streamsize const sSize = std::cout.precision();
          std::stringstream ss;
          ss << std::setprecision(12) << interpolationVector[nodeBirth] << " > "
             << interpolationVector[nodeDeath];
          printErr(ss.str());
          std::cout.precision(sSize);
        }
        // - Verify nesting condition
        bool const verifyNesting = !baryTree->notNeedToNormalize(i);
        if(normalizedWasserstein_ and verifyNesting) {
          if(interpolationVector[nodeBirth] < interParentBirth
             and (interParentBirth - interpolationVector[nodeBirth]) < eps)
            interpolationVector[nodeBirth] = interParentBirth;
          if(interpolationVector[nodeBirth] > interParentDeath
             and (interpolationVector[nodeBirth] - interParentDeath) < eps) {
            interpolationVector[nodeBirth] = interParentDeath;
          }
          if(interpolationVector[nodeDeath] < interParentBirth
             and (interParentBirth - interpolationVector[nodeDeath]) < eps) {
            interpolationVector[nodeDeath] = interParentBirth;
          }
          if(interpolationVector[nodeDeath] > interParentDeath
             and (interpolationVector[nodeDeath] - interParentDeath) < eps)
            interpolationVector[nodeDeath] = interParentDeath;
        }
        if(normalizedWasserstein_ and verifyNesting
           and (interpolationVector[nodeBirth] < interParentBirth
                or interpolationVector[nodeBirth] > interParentDeath
                or interpolationVector[nodeDeath] < interParentBirth
                or interpolationVector[nodeDeath] > interParentDeath)) {
          ++cptNotNesting;
          printMsg(ssT.str());
          std::streamsize const sSize = std::cout.precision();
          std::stringstream ss;
          ss << std::setprecision(12) << interpolationVector[nodeBirth] << " _ "
             << interpolationVector[nodeDeath] << " --- " << interParentBirth
             << " _ " << interParentDeath;
          printErr(ss.str());
          std::cout.precision(sSize);
        }

        // - Push children to the queue
        std::vector<ftm::idNode> children;
        baryTree->getChildren(i, children);
        for(auto child : children)
          queue.emplace(child);
      }
      if(noNestingShortener != 0)
        printWrn("[getInterpolationVector] "
                 + std::to_string(noNestingShortener)
                 + " nesting vector shortener.");
      if(noDiagonalShortener != 0)
        printMsg("[getInterpolationVector] "
                   + std::to_string(noDiagonalShortener)
                   + " diagonal vector shortener.",
                 debug::Priority::DETAIL);
      if(cptBelowDiagonal != 0)
        printErr("[getInterpolationVector] point below diagonal.");
      if(cptNotNesting != 0)
        printErr("[getInterpolationVector] nesting condition not ok.");

      interpolationVectorT.resize(scalarsVector.size());
      for(unsigned int i = 0; i < interpolationVectorT.size(); ++i)
        interpolationVectorT[i] = interpolationVector[i];
    }

    template <class dataType>
    void getInterpolation(ftm::MergeTree<dataType> &barycenter,
                          std::vector<double *> &v,
                          std::vector<double *> &v2,
                          size_t vSize,
                          double t,
                          ftm::MergeTree<dataType> &interpolated,
                          bool transposeVector = true) {
      ftm::FTMTree_MT *barycenterTree = &(barycenter.tree);

      // Compute new scalars vector
      std::vector<dataType> newScalarsVector;
      getInterpolationVector<dataType>(
        barycenter, v, v2, vSize, t, newScalarsVector, transposeVector);

      // Create new tree
      interpolated
        = ttk::ftm::createEmptyMergeTree<dataType>(newScalarsVector.size());
      ttk::ftm::setTreeScalars<dataType>(interpolated, newScalarsVector);
      ftm::FTMTree_MT *treeNew = &(interpolated.tree);

      // Copy the barycenter tree structure
      treeNew->copyMergeTreeStructure(barycenterTree);

      // Delete diagonal and almost-diagonal points
      persistenceThresholding<dataType>(treeNew, 0.001); // 0.001 %
    }

    template <class dataType>
    void getInterpolation(ftm::MergeTree<dataType> &barycenter,
                          std::vector<std::vector<double>> &v,
                          std::vector<std::vector<double>> &v2,
                          double t,
                          ftm::MergeTree<dataType> &interpolated) {
      std::vector<double *> pV, pV2;
      vectorsToPointers(v, pV);
      vectorsToPointers(v2, pV2);
      return getInterpolation<dataType>(
        barycenter, pV, pV2, v[0].size(), t, interpolated, false);
    }

    template <class dataType>
    void getMultiInterpolation(ftm::MergeTree<dataType> &barycenter,
                               std::vector<std::vector<double *>> &vS,
                               std::vector<std::vector<double *>> &v2s,
                               size_t vSize,
                               std::vector<double> &ts,
                               ftm::MergeTree<dataType> &interpolated,
                               bool transposeVector = true) {
      getInterpolation<dataType>(
        barycenter, vS[0], v2s[0], vSize, ts[0], interpolated, transposeVector);
      for(unsigned int i = 1; i < vS.size(); ++i) {
        ftm::MergeTree<dataType> newInterpolated;
        getInterpolation<dataType>(interpolated, vS[i], v2s[i], vSize, ts[i],
                                   newInterpolated, transposeVector);
        interpolated = newInterpolated;
      }
    }

    template <class dataType>
    void
      getMultiInterpolation(ftm::MergeTree<dataType> &barycenter,
                            std::vector<std::vector<std::vector<double>>> &vS,
                            std::vector<std::vector<std::vector<double>>> &v2s,
                            std::vector<double> &ts,
                            ftm::MergeTree<dataType> &interpolated) {
      std::vector<std::vector<double *>> pVS, pV2s;
      vectorOfVectorsToPointers(vS, pVS);
      vectorOfVectorsToPointers(v2s, pV2s);
      getMultiInterpolation(
        barycenter, pVS, pV2s, vS[0][0].size(), ts, interpolated, false);
    }

    template <class dataType>
    void
      getMultiInterpolation(ftm::MergeTree<dataType> &barycenter,
                            std::vector<std::vector<std::vector<double>>> &vS,
                            std::vector<std::vector<std::vector<double>>> &v2s,
                            std::vector<std::vector<double>> &v,
                            std::vector<std::vector<double>> &v2,
                            std::vector<double> &ts,
                            double t,
                            ftm::MergeTree<dataType> &interpolated) {
      std::vector<std::vector<std::vector<double>>> vSTemp = vS, v2sTemp = v2s;
      vSTemp.push_back(v);
      v2sTemp.push_back(v2);
      std::vector<double> tsTemp = ts;
      if(tsTemp.size() != 0)
        tsTemp[vSTemp.size() - 1] = t;
      else
        tsTemp.push_back(t);
      getMultiInterpolation(barycenter, vSTemp, v2sTemp, tsTemp, interpolated);
    }

    // ----------------------------------------------------------------------------
    // Vector Utils
    // ----------------------------------------------------------------------------
    void callGramSchmidt(std::vector<std::vector<double>> &vS,
                         std::vector<double> &v,
                         std::vector<double> &newV);

    void vectorToPointer(std::vector<double> &vec, double *&pVec);

    void vectorsToPointers(std::vector<std::vector<double>> &vec,
                           std::vector<double *> &pVec);

    void vectorOfVectorsToPointers(
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<double *>> &pVS);

    void pointerToVector(double *pVec, size_t size, std::vector<double> &vec);

    void pointersToVectors(std::vector<double *> &pVec,
                           std::vector<size_t> sizes,
                           std::vector<std::vector<double>> &vec);

    void pointersToVectors(std::vector<double *> &pVec,
                           size_t size,
                           std::vector<std::vector<double>> &vec);
  }; // MergeTreePrincipalGeodesicsBase class

} // namespace ttk
