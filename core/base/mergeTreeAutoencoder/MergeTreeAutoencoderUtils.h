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

    /**
     * @brief Fix the scalars of a merge tree to ensure that the nesting
     * condition is respected.
     *
     * @param[in] mTree Merge tree to process.
     */
    void fixTreePrecisionScalars(ftm::MergeTree<float> &mTree);

    /**
     * @brief Fix the scalars of a merge tree to ensure that the nesting
     * condition is respected.
     *
     * @param[in] scalarsVector scalars array to process.
     * @param[in] node node to adjust.
     * @param[in] refNode reference node.
     */
    void adjustNestingScalars(std::vector<float> &scalarsVector,
                              ftm::idNode node,
                              ftm::idNode refNode);

    /**
     * @brief Create a balanced BDT structure (for output basis initialization).
     *
     * @param[in] parents vector containing the possible parents for each node.
     * @param[in] children vector containing the possible children for each
     * node.
     * @param[in] scalarsVector vector containing the scalars value.
     * @param[out] childrenFinal output vector containing the children of each
     * node, representing the tree structure.
     * @param[in] threadNumber number of threads for parallel sort.
     */
    void createBalancedBDT(std::vector<std::vector<ftm::idNode>> &parents,
                           std::vector<std::vector<ftm::idNode>> &children,
                           std::vector<float> &scalarsVector,
                           std::vector<std::vector<ftm::idNode>> &childrenFinal,
                           int threadNumber = 1);

    /**
     * @brief Util function to print pairs of a merge tree.
     *
     * @param[in] mTree merge tree to process.
     * @param[in] useBD if the merge tree is in branch decomposition mode or
     * not.
     */
    void printPairs(ftm::MergeTree<float> &mTree, bool useBD = true);

  } // namespace wae

} // namespace ttk
