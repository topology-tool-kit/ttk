/// \ingroup base
/// \class ttk::FTMTreePPUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEPPUTILS_H
#define _FTMTREEPPUTILS_H

#pragma once

#include <FTMTree.h>
#include <FTMTreePP.h>
#include <FTMTreeUtils.h>

#include <ttkUtils.h>

using namespace ttk;
using namespace ftm;

template <class dataType>
void getPersistencePairs(
  ftm::FTMTree_MT *tree,
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &pairs) {
  ttk::ftm::FTMTreePP pairsCompute;
  pairsCompute.setCustomTree(tree);
  pairsCompute.computePersistencePairs<dataType>(
    pairs, isJoinTree<dataType>(tree));
}

template <class dataType>
std::vector<std::tuple<SimplexId, SimplexId, dataType>>
  computePersistencePairs(ftm::FTMTree_MT *tree) {
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> pairs;
  getPersistencePairs<dataType>(tree, pairs);
  for(auto pair : pairs) {
    if(tree->getNode(std::get<0>(pair))->getOrigin() < std::get<0>(pair)
       and tree->getNode(std::get<0>(pair))->getOrigin() >= 0)
      tree->getNode(tree->getNode(std::get<0>(pair))->getOrigin())
        ->setOrigin(std::get<1>(pair));

    tree->getNode(std::get<0>(pair))->setOrigin(std::get<1>(pair));
    tree->getNode(std::get<1>(pair))->setOrigin(std::get<0>(pair));
  }
  return pairs;
}

#endif
