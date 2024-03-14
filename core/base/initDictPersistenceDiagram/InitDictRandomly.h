#pragma once

#include <PersistenceDiagramUtils.h>
#include <Wrapper.h>

#include <algorithm>
#include <array>
#include <tuple>

namespace ttk {
  using Matrix = std::vector<std::vector<double>>;

  class InitRandomDict : public Debug {

  public:
    InitRandomDict() {
      this->setDebugMsgPrefix("InitRandomDict");
    };

    void execute(std::vector<ttk::DiagramType> &DictDiagrams,
                 const std::vector<ttk::DiagramType> &datas,
                 const int nbAtoms,
                 const int seed);
  };
} // namespace ttk
