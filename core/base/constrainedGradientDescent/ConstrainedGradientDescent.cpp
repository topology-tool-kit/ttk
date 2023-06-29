#include <ConstrainedGradientDescent.h>
#include <cmath>
#include <csignal>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif // TTK_ENABLE_EIGEN

using namespace ttk;

void ConstrainedGradientDescent::executeWeightsProjected(
  std::vector<Matrix> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  bool maxEigenValue) {
  gradientDescentWeights(hessianList, weights, grad, maxEigenValue);
  projectionOnSimplex(weights);
}

void ConstrainedGradientDescent::executeAtoms(
  std::vector<ttk::DiagramType> &DictDiagrams,
  const std::vector<std::vector<ttk::MatchingType>> &matchings,
  const ttk::DiagramType &Barycenter,
  const std::vector<Matrix> &gradsLists,
  const std::vector<int> &checkerAtomsExt,
  std::vector<std::vector<int>> &projForDiag,
  ttk::DiagramType &featuresToAdd,
  std::vector<std::array<double, 2>> &projLocations,
  std::vector<std::vector<double>> &vectorForProjContrib,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  ttk::DiagramType &infoToAdd) {
  gradientDescentAtoms(DictDiagrams, matchings, Barycenter, gradsLists,
                       checkerAtomsExt, projForDiag, featuresToAdd,
                       projLocations, vectorForProjContrib, pairToAddGradList,
                       infoToAdd);
}

// simple projection on simplex, aka where a vector has positive elements and
// sum to 1.
void ConstrainedGradientDescent::projectionOnSimplex(
  std::vector<double> &weights) {

  int n = weights.size();
  std::vector<double> copy_temp = weights;
  std::sort(copy_temp.rbegin(), copy_temp.rend());
  double K = 1.;
  double somme_u = copy_temp[0];
  double theta = (somme_u - 1.) / K;
  while(K < n && (somme_u + copy_temp[K] - 1.) / (K + 1.) < copy_temp[K]) {
    somme_u += copy_temp[K];
    K += 1.;
    theta = (somme_u - 1.) / K;
  }
  for(int i = 0; i < n; ++i) {
    weights[i] = std::max(weights[i] - theta, 0.);
  }

  double sum = 0.;
  for(int i = 0; i < n - 1; ++i) {
    weights[i] = trunc(weights[i] * 1e8) / 1e8;
    sum += weights[i];
  }
  weights[n - 1] = 1. - sum;
}

void ConstrainedGradientDescent::gradientDescentWeights(
  std::vector<Matrix> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  bool maxEigenValue) {

  int n = weights.size();
  double stepWeight;
  double L = 0.;
#ifndef TTK_ENABLE_EIGEN
  maxEigenValue = false;
#endif // TTK_ENABLE_EIGEN
  if(maxEigenValue) {
#ifdef TTK_ENABLE_EIGEN
    for(size_t i = 0; i < hessianList.size(); ++i) {
      auto &hessian = hessianList[i];
      int m = hessian.size();
      Eigen::MatrixXd H(m, m);
      for(size_t j = 0; j < hessian.size(); ++j) {
        for(size_t k = 0; k < hessian.size(); ++k) {
          H(j, k) = hessian[j][k];
        }
      }
      Eigen::EigenSolver<Eigen::MatrixXd> es;
      es.compute(H, false);
      Eigen::VectorXcd eigvals = es.eigenvalues();
      L += eigvals.lpNorm<Eigen::Infinity>();
    }
#endif // TTK_ENABLE_EIGEN
  } else {
    for(size_t i = 0; i < hessianList.size(); ++i) {
      auto &hessian = hessianList[i];
      for(size_t k = 0; k < hessian.size(); ++k) {
        double diag = hessian[k][k];
        L += 1. * diag;
      }
    }
  }

  if(L > 0) {
    stepWeight = 1. / L;
  } else {
    stepWeight = 0;
  }

  for(int i = 0; i < n; ++i) {
    weights[i] = weights[i] - stepWeight * grad[i];
  }
}

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// METTRE TIMER POUR VOIR QUOI PARALELLISER
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

void ConstrainedGradientDescent::gradientDescentAtoms(
  std::vector<ttk::DiagramType> &DictDiagrams,
  const std::vector<std::vector<ttk::MatchingType>> &matchings,
  const ttk::DiagramType &Barycenter,
  const std::vector<Matrix> &gradsLists,
  const std::vector<int> &checkerAtomsExt,
  std::vector<std::vector<int>> &projForDiag,
  ttk::DiagramType &featuresToAdd,
  std::vector<std::array<double, 2>> &projLocations,
  std::vector<std::vector<double>> &vectorForProjContrib,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  ttk::DiagramType &infoToAdd) {

  // Here vector of diagramTuple because it is not a persistence diagram per
  // say.
  // we get the right pairs to update for each barycenter pair.

  std::vector<double> miniBirth(matchings.size());
  for(size_t i = 0; i < matchings.size(); ++i) {
    auto &t = DictDiagrams[i][0];
    miniBirth[i] = t.birth.sfValue;
  }

  std::vector<std::vector<std::array<double, 2>>> gradBuffersList(
    Barycenter.size());
  std::vector<std::vector<double>> projectionsBuffer(Barycenter.size());

  for(size_t i = 0; i < Barycenter.size(); ++i) {
    projectionsBuffer[i].resize(matchings.size());
  }

  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    gradBuffersList[i].resize(matchings.size());
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  for(size_t i = 0; i < checker.size(); ++i) {
    checker[i].resize(matchings.size());
  }
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<std::vector<int>> trackerDiagonal(Barycenter.size());
  std::vector<std::vector<int>> trackerMatch(Barycenter.size());
  for(size_t i = 0; i < trackerDiagonal.size(); ++i) {
    trackerDiagonal[i].resize(matchings.size());
  }
  for(size_t i = 0; i < trackerMatch.size(); ++i) {
    trackerMatch[i].resize(matchings.size());
  }

  for(size_t i = 0; i < matchings.size(); ++i) {
    for(size_t j = 0; j < matchings[i].size(); ++j) {
      const auto &t = matchings[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      if(Id2 < 0 || static_cast<SimplexId>(gradBuffersList.size()) <= Id2
         || static_cast<SimplexId>(DictDiagrams[i].size()) <= Id1) {
        continue;
      } else {
        if(Id1 < 0) {
          const auto &t3 = Barycenter[Id2];
          auto &point = gradBuffersList[Id2][i];
          const double birthBarycenter = t3.birth.sfValue;
          const double deathBarycenter = t3.death.sfValue;
          const double birthDeathAtom
            = birthBarycenter + (deathBarycenter - birthBarycenter) / 2.;
          point[0] = birthDeathAtom;
          point[1] = birthDeathAtom;
          checker[Id2][i] = i;
          if(checker[Id2].size() > matchings.size()) {
            std::raise(SIGINT);
          }
          tracker[Id2] = 1;
          trackerMatch[Id2][i] = Id1;
          if(static_cast<SimplexId>(DictDiagrams[i].size()) <= Id1) {
            std::cout << "ID1: " << Id1 << std::endl;
          }
          trackerDiagonal[Id2][i] = 1;
          projectionsBuffer[Id2][i] = birthDeathAtom;

        } else {
          const auto &t2 = DictDiagrams[i][Id1];
          auto &point = gradBuffersList[Id2][i];
          const double birthAtom = t2.birth.sfValue;
          const double deathAtom = t2.death.sfValue;
          point[0] = birthAtom;
          point[1] = deathAtom;
          checker[Id2][i] = i;
          if(checker[Id2].size() > matchings.size()) {
            std::raise(SIGINT);
          }
          tracker[Id2] = 1;
          trackerMatch[Id2][i] = Id1;
          trackerDiagonal[Id2][i] = 0;
          projectionsBuffer[Id2][i] = 0.;
        }
      }
    }
  }

  gradBuffersList.insert(
    gradBuffersList.end(), pairToAddGradList.begin(), pairToAddGradList.end());
  for(size_t j = 0; j < pairToAddGradList.size(); ++j) {
    std::vector<int> temp1(matchings.size(), 1);
    std::vector<int> temp2(matchings.size(), -1);
    std::vector<int> temp3(matchings.size());
    for(size_t l = 0; l < matchings.size(); ++l) {
      temp3[l] = static_cast<int>(l);
    }
    trackerDiagonal.push_back(temp1);
    trackerMatch.push_back(temp1);
    checker.push_back(temp3);
    tracker.push_back(1);
  }

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {

      continue;
    } else {

      std::vector<double> pos(gradBuffersList[i].size(), 0.);
      for(size_t j = 0; j < checker[i].size(); ++j) {
        auto &t = gradBuffersList[i][checker[i][j]];
        double birth = t[0];
        double death = t[1];
        pos[j] = death - birth;
      }
      std::vector<bool> pos2(pos.size(), false);
      std::vector<double> temp2;
      for(size_t p = 0; p < checker[i].size(); ++p) {
        auto &t = gradBuffersList[i][checker[i][p]];
        double birth = t[0];
        pos2[p] = birth == 0.;
        if(birth > 0.) {
          temp2.push_back(birth);
        }
      }
      for(size_t p = 0; p < checker[i].size(); ++p) {
        if(pos2[p]) {

          auto &t = gradBuffersList[i][checker[i][p]];

          t[1] = t[1] - (this->stepAtom) * gradsLists[i][checker[i][p]][1];
        } else if(pos[p] < 1e-7) {
          continue;
        } else {

          auto &t0 = gradBuffersList[i][checker[i][p]];

          t0[0] = t0[0] - (this->stepAtom) * gradsLists[i][checker[i][p]][0];
          t0[1] = t0[1] - (this->stepAtom) * gradsLists[i][checker[i][p]][1];

          if(t0[0] > t0[1]) {
            t0[1] = t0[0];
          }
        }
      }
    }
  }

  for(size_t i = 0; i < checker.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {
      continue;
    } else {
      if(i < Barycenter.size()) {
        for(size_t j = 0; j < checker[i].size(); ++j) {
          auto &trackerTemp = trackerMatch[i][j];
          if(trackerDiagonal[i][j] == 1 || trackerTemp == -1
             || trackerTemp > 10000) {
            std::vector<int> projAndIndex(DictDiagrams.size() + 1);
            for(size_t m = 0; m < DictDiagrams.size(); ++m) {
              projAndIndex[m] = trackerMatch[i][m];
            }

          } else {
            auto &index = checker[i][j];
            auto &t2 = gradBuffersList[i][index];

            auto &t1 = DictDiagrams[index][trackerTemp];

            if(t2[0] < miniBirth[index]) {
              t1.birth.sfValue = miniBirth[index];
            } else {
              t1.birth.sfValue = t2[0];
            }
            if(t2[1] < t2[0]) {
              t1.birth.sfValue = t2[0];
              t1.death.sfValue = t2[0];
            } else {
              t1.death.sfValue = t2[1];
            }
          }
        }
      } else {
        const auto &infos = infoToAdd[static_cast<int>(i)
                                      - static_cast<int>(Barycenter.size())];

        for(size_t j = 0; j < checker[i].size(); ++j) {
          auto &trackerTemp = trackerMatch[i][j];
          if(trackerDiagonal[i][j] == 1 || trackerTemp == -1
             || trackerTemp > 10000) {
            auto &index = checker[i][j];
            std::vector<int> projAndIndex(DictDiagrams.size() + 1);
            int atomIndex = static_cast<int>(index);
            for(size_t m = 0; m < DictDiagrams.size(); ++m) {
              projAndIndex[m] = trackerMatch[i][m];
            }
            projAndIndex[DictDiagrams.size()] = atomIndex;
            projForDiag.push_back(projAndIndex);
            featuresToAdd.push_back(infos);
            projLocations.push_back(gradBuffersList[i][index]);
            vectorForProjContrib.push_back(gradsLists[i][index]);
          } else {
            auto &index = checker[i][j];
            auto &t2 = gradBuffersList[i][index];
            auto &t1 = DictDiagrams[index][trackerTemp];
            if(t2[0] < miniBirth[index]) {
              t1.birth.sfValue = miniBirth[index];
            } else {
              t1.birth.sfValue = t2[0];
            }
            if(t2[1] < t2[0]) {
              t1.birth.sfValue = t2[0];
              t1.death.sfValue = t2[0];
            } else {
              t1.death.sfValue = t2[1];
            }
          }
        }
      }
    }
  }
}

void ConstrainedGradientDescent::setStep(double factEquiv) {
  this->stepAtom = 1. / (2. * 2. * factEquiv);
}

void ConstrainedGradientDescent::reduceStep() {
  this->stepAtom = this->stepAtom / 2.;
}