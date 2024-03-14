#pragma once

// #include <PDClustering.h>
// #include <PersistenceDiagramBarycenter.h>
//  #include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>
#include <Wrapper.h>

#include <algorithm>
#include <array>

namespace ttk {
  using Matrix = std::vector<std::vector<double>>;

  class InitDictPersistenceDiagram : public Debug {

  public:
    InitDictPersistenceDiagram() {
      this->setDebugMsgPrefix("InitFarBorderDict");
    };

    void execute(std::vector<ttk::DiagramType> &DictDiagrams,
                 const std::vector<ttk::DiagramType> &datas,
                 const int nbAtoms,
                 bool do_min_,
                 bool do_sad_,
                 bool do_max_);

  protected:
    int getNextIndex(const Matrix &distMatrix,
                     const std::vector<int> &indices) const;

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    double Lambda{0};
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
  };
} // namespace ttk
