/// \ingroup base
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// \brief Provide utils methods related to Merge Trees and Persistence Diagrams
/// Auto-Encoders.

#pragma once

#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <FTMTree_MT.h>

namespace ttk {

  namespace wae {

    void fixTreePrecisionScalars(ftm::MergeTree<float> &mTree);

    void adjustNestingScalars(std::vector<float> &scalarsVector,
                              ftm::idNode node,
                              ftm::idNode refNode);

    void
      createBalancedBDT(std::vector<std::vector<ftm::idNode>> &parents,
                        std::vector<std::vector<ftm::idNode>> &children,
                        std::vector<float> &scalarsVector,
                        std::vector<std::vector<ftm::idNode>> &childrenFinal);

    void printPairs(ftm::MergeTree<float> &mTree, bool useBD = true);

  } // namespace wae

} // namespace ttk
