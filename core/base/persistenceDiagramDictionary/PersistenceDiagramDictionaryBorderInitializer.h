#pragma once

#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>
#include <Wrapper.h>

namespace ttk {
  using Matrix = std::vector<std::vector<double>>;

  class PersistenceDiagramDictionaryBorderInitializer : public Debug {

  public:
    PersistenceDiagramDictionaryBorderInitializer() {
      this->setDebugMsgPrefix("InitFarBorderDict");
    };

    void execute(std::vector<ttk::DiagramType> &DictDiagrams,
                 const std::vector<ttk::DiagramType> &datas,
                 const int &nbAtoms,
                 bool do_min_,
                 bool do_sad_,
                 bool do_max_);

  protected:
    int getNextIndex(const Matrix &distMatrix,
                     const std::vector<int> &indices) const;
  };
} // namespace ttk
