#pragma once

#include <PersistenceDiagramUtils.h>
#include <Wrapper.h>

#include <algorithm>
#include <array>
#include <tuple>

namespace ttk {
  using Matrix = std::vector<std::vector<double>>;

  class ConstrainedGradientDescent : public Debug {

  public:
    ConstrainedGradientDescent() {
      this->setDebugMsgPrefix("ConstrainedGradientDescent");
    };

    void executeWeightsProjected(std::vector<Matrix> &hessianList,
                                 std::vector<double> &weights,
                                 const std::vector<double> &grad,
                                 bool maxEigenValue);

    void executeAtoms(
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
      ttk::DiagramType &infoToAdd);

    void setStep(double factEquiv);
    void reduceStep();
    // void executeAtoms(std::vector<Diagram> &DictDiagrams);

    // inline void setNbAtoms(const int nbAtoms) {
    // NbAtoms = nbAtoms;
    //}

  protected:
    void projectionOnSimplex(std::vector<double> &weights);

    void gradientDescentWeights(std::vector<Matrix> &hessianList,
                                std::vector<double> &weights,
                                const std::vector<double> &grad,
                                bool maxEigenValue);

    void gradientDescentAtoms(
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
      ttk::DiagramType &infoToAdd);

    double stepAtom;
  };

} // namespace ttk
