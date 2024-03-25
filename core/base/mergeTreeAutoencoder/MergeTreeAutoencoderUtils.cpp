#include <MergeTreeAutoencoderUtils.h>

void ttk::wae::fixTreePrecisionScalars(ftm::MergeTree<float> &mTree) {
  double eps = 1e-6;
  auto shiftSubtree
    = [&mTree, &eps](ftm::idNode node, ftm::idNode birthNodeParent,
                     ftm::idNode deathNodeParent, std::vector<float> &scalars,
                     bool invalidBirth, bool invalidDeath) {
        std::queue<ftm::idNode> queue;
        queue.emplace(node);
        while(!queue.empty()) {
          ftm::idNode nodeT = queue.front();
          queue.pop();
          auto birthDeathNode = mTree.tree.getBirthDeathNode<float>(node);
          auto birthNode = std::get<0>(birthDeathNode);
          auto deathNode = std::get<1>(birthDeathNode);
          if(invalidBirth)
            scalars[birthNode] = scalars[birthNodeParent] + 2 * eps;
          if(invalidDeath)
            scalars[deathNode] = scalars[deathNodeParent] - 2 * eps;
          std::vector<ftm::idNode> children;
          mTree.tree.getChildren(nodeT, children);
          for(auto &child : children)
            queue.emplace(child);
        }
      };
  std::vector<float> scalars;
  getTreeScalars(mTree, scalars);
  std::queue<ftm::idNode> queue;
  auto root = mTree.tree.getRoot();
  queue.emplace(root);
  while(!queue.empty()) {
    ftm::idNode node = queue.front();
    queue.pop();
    auto birthDeathNode = mTree.tree.getBirthDeathNode<float>(node);
    auto birthNode = std::get<0>(birthDeathNode);
    auto deathNode = std::get<1>(birthDeathNode);
    auto birthDeathNodeParent
      = mTree.tree.getBirthDeathNode<float>(mTree.tree.getParentSafe(node));
    auto birthNodeParent = std::get<0>(birthDeathNodeParent);
    auto deathNodeParent = std::get<1>(birthDeathNodeParent);
    bool invalidBirth = (scalars[birthNode] <= scalars[birthNodeParent] + eps);
    bool invalidDeath = (scalars[deathNode] >= scalars[deathNodeParent] - eps);
    if(!mTree.tree.isRoot(node) and (invalidBirth or invalidDeath))
      shiftSubtree(node, birthNodeParent, deathNodeParent, scalars,
                   invalidBirth, invalidDeath);
    std::vector<ftm::idNode> children;
    mTree.tree.getChildren(node, children);
    for(auto &child : children)
      queue.emplace(child);
  }
  ftm::setTreeScalars<float>(mTree, scalars);
}

void ttk::wae::adjustNestingScalars(std::vector<float> &scalarsVector,
                                    ftm::idNode node,
                                    ftm::idNode refNode) {
  float birth = scalarsVector[refNode * 2];
  float death = scalarsVector[refNode * 2 + 1];
  auto getSign = [](float v) { return (v > 0 ? 1 : -1); };
  auto getPrecValue = [&getSign](float v, bool opp = false) {
    return v * (1 + (opp ? -1 : 1) * getSign(v) * 1e-6);
  };
  // Shift scalars
  if(scalarsVector[node * 2 + 1] > getPrecValue(death, true)) {
    float diff = scalarsVector[node * 2 + 1] - getPrecValue(death, true);
    scalarsVector[node * 2] -= diff;
    scalarsVector[node * 2 + 1] -= diff;
  } else if(scalarsVector[node * 2] < getPrecValue(birth)) {
    float diff = getPrecValue(birth) - scalarsVector[node * 2];
    scalarsVector[node * 2] += getPrecValue(diff);
    scalarsVector[node * 2 + 1] += getPrecValue(diff);
  }
  // Cut scalars
  if(scalarsVector[node * 2] < getPrecValue(birth))
    scalarsVector[node * 2] = getPrecValue(birth);
  if(scalarsVector[node * 2 + 1] > getPrecValue(death, true))
    scalarsVector[node * 2 + 1] = getPrecValue(death, true);
}

void ttk::wae::createBalancedBDT(
  std::vector<std::vector<ftm::idNode>> &parents,
  std::vector<std::vector<ftm::idNode>> &children,
  std::vector<float> &scalarsVector,
  std::vector<std::vector<ftm::idNode>> &childrenFinal,
  int threadNumber) {
  // ----- Some variables
  unsigned int noNodes = scalarsVector.size() / 2;
  childrenFinal.resize(noNodes);
  int mtLevel = ceil(log(noNodes * 2) / log(2)) + 1;
  int bdtLevel = mtLevel - 1;
  int noDim = bdtLevel;

  // ----- Get node levels
  std::vector<int> nodeLevels(noNodes, -1);
  std::queue<ftm::idNode> queueLevels;
  std::vector<int> noChildDone(noNodes, 0);
  for(unsigned int i = 0; i < children.size(); ++i) {
    if(children[i].size() == 0) {
      queueLevels.emplace(i);
      nodeLevels[i] = 1;
    }
  }
  while(!queueLevels.empty()) {
    ftm::idNode node = queueLevels.front();
    queueLevels.pop();
    for(auto &parent : parents[node]) {
      ++noChildDone[parent];
      nodeLevels[parent] = std::max(nodeLevels[parent], nodeLevels[node] + 1);
      if(noChildDone[parent] >= (int)children[parent].size())
        queueLevels.emplace(parent);
    }
  }

  // ----- Sort heuristic lambda
  auto sortChildren = [&parents, &scalarsVector, &noNodes, &threadNumber](
                        ftm::idNode nodeOrigin, std::vector<bool> &nodeDone,
                        std::vector<std::vector<ftm::idNode>> &childrenT) {
    double refPers = scalarsVector[1] - scalarsVector[0];
    auto getRemaining = [&nodeDone](std::vector<ftm::idNode> &vec) {
      unsigned int remaining = 0;
      for(auto &e : vec)
        remaining += (not nodeDone[e]);
      return remaining;
    };
    std::vector<unsigned int> parentsRemaining(noNodes, 0),
      childrenRemaining(noNodes, 0);
    for(auto &child : childrenT[nodeOrigin]) {
      parentsRemaining[child] = getRemaining(parents[child]);
      childrenRemaining[child] = getRemaining(childrenT[child]);
    }
    TTK_FORCE_USE(threadNumber);
    TTK_PSORT(
      threadNumber, childrenT[nodeOrigin].begin(), childrenT[nodeOrigin].end(),
      [&](ftm::idNode nodeI, ftm::idNode nodeJ) {
        double persI = scalarsVector[nodeI * 2 + 1] - scalarsVector[nodeI * 2];
        double persJ = scalarsVector[nodeJ * 2 + 1] - scalarsVector[nodeJ * 2];
        return parentsRemaining[nodeI] + childrenRemaining[nodeI]
                 - persI / refPers * noNodes
               < parentsRemaining[nodeJ] + childrenRemaining[nodeJ]
                   - persJ / refPers * noNodes;
      });
  };

  // ----- Greedy approach to find balanced BDT structures
  const auto findStructGivenDim =
    [&children, &noNodes, &nodeLevels](
      ftm::idNode _nodeOrigin, int _dimToFound, bool _searchMaxDim,
      std::vector<bool> &_nodeDone, std::vector<bool> &_dimFound,
      std::vector<std::vector<ftm::idNode>> &_childrenFinalOut) {
      // --- Recursive lambda
      auto findStructGivenDimImpl =
        [&children, &noNodes, &nodeLevels](
          ftm::idNode nodeOrigin, int dimToFound, bool searchMaxDim,
          std::vector<bool> &nodeDone, std::vector<bool> &dimFound,
          std::vector<std::vector<ftm::idNode>> &childrenFinalOut,
          auto &findStructGivenDimRef) mutable {
          childrenFinalOut.resize(noNodes);
          // - Find structures
          int dim = (searchMaxDim ? dimToFound - 1 : 0);
          unsigned int i = 0;
          //
          auto searchMaxDimReset = [&i, &dim, &nodeDone]() {
            --dim;
            i = 0;
            unsigned int noDone = 0;
            for(auto done : nodeDone)
              if(done)
                ++noDone;
            return noDone == nodeDone.size() - 1; // -1 for root
          };
          while(i < children[nodeOrigin].size()) {
            auto child = children[nodeOrigin][i];
            // Skip if child was already processed
            if(nodeDone[child]) {
              // If we have processed all children while searching for max
              // dim then restart at the beginning to find a lower dim
              if(searchMaxDim and i == children[nodeOrigin].size() - 1) {
                if(searchMaxDimReset())
                  break;
              } else
                ++i;
              continue;
            }
            if(dim == 0) {
              // Base case
              childrenFinalOut[nodeOrigin].emplace_back(child);
              nodeDone[child] = true;
              dimFound[0] = true;
              if(dimToFound <= 1 or searchMaxDim)
                return true;
              ++dim;
            } else {
              // General case
              std::vector<std::vector<ftm::idNode>> childrenFinalDim;
              std::vector<bool> nodeDoneDim;
              std::vector<bool> dimFoundDim(dim);
              bool found = false;
              if(nodeLevels[child] > dim) {
                nodeDoneDim = nodeDone;
                found = findStructGivenDimRef(child, dim, false, nodeDoneDim,
                                              dimFoundDim, childrenFinalDim,
                                              findStructGivenDimRef);
              }
              if(found) {
                dimFound[dim] = true;
                childrenFinalOut[nodeOrigin].emplace_back(child);
                for(unsigned int j = 0; j < childrenFinalDim.size(); ++j)
                  for(auto &e : childrenFinalDim[j])
                    childrenFinalOut[j].emplace_back(e);
                nodeDone[child] = true;
                for(unsigned int j = 0; j < nodeDoneDim.size(); ++j)
                  nodeDone[j] = nodeDone[j] || nodeDoneDim[j];
                // Return if it is the last dim to found
                if(dim == dimToFound - 1 and not searchMaxDim)
                  return true;
                // Reset index if we search for the maximum dim
                if(searchMaxDim) {
                  if(searchMaxDimReset())
                    break;
                } else {
                  ++dim;
                }
                continue;
              } else if(searchMaxDim and i == children[nodeOrigin].size() - 1) {
                // If we have processed all children while searching for max
                // dim then restart at the beginning to find a lower dim
                if(searchMaxDimReset())
                  break;
                continue;
              }
            }
            ++i;
          }
          return false;
        };
      return findStructGivenDimImpl(_nodeOrigin, _dimToFound, _searchMaxDim,
                                    _nodeDone, _dimFound, _childrenFinalOut,
                                    findStructGivenDimImpl);
    };
  std::vector<bool> dimFound(noDim - 1, false);
  std::vector<bool> nodeDone(noNodes, false);
  for(unsigned int i = 0; i < children.size(); ++i)
    sortChildren(i, nodeDone, children);
  Timer t_find;
  ftm::idNode startNode = 0;
  findStructGivenDim(startNode, noDim, true, nodeDone, dimFound, childrenFinal);

  // ----- Greedy approach to create non found structures
  const auto createStructGivenDim =
    [&children, &noNodes, &findStructGivenDim, &nodeLevels](
      int _nodeOrigin, int _dimToCreate, std::vector<bool> &_nodeDone,
      ftm::idNode &_structOrigin, std::vector<float> &_scalarsVectorOut,
      std::vector<std::vector<ftm::idNode>> &_childrenFinalOut) {
      // --- Recursive lambda
      auto createStructGivenDimImpl =
        [&children, &noNodes, &findStructGivenDim, &nodeLevels](
          int nodeOrigin, int dimToCreate, std::vector<bool> &nodeDoneImpl,
          ftm::idNode &structOrigin, std::vector<float> &scalarsVectorOut,
          std::vector<std::vector<ftm::idNode>> &childrenFinalOut,
          auto &createStructGivenDimRef) mutable {
          // Deduction of auto lambda type
          if(false)
            return;
          // - Find structures of lower dimension
          int dimToFound = dimToCreate - 1;
          std::vector<std::vector<std::vector<ftm::idNode>>> childrenFinalT(2);
          std::array<ftm::idNode, 2> structOrigins;
          for(unsigned int n = 0; n < 2; ++n) {
            bool found = false;
            for(unsigned int i = 0; i < children[nodeOrigin].size(); ++i) {
              auto child = children[nodeOrigin][i];
              if(nodeDoneImpl[child])
                continue;
              if(dimToFound != 0) {
                if(nodeLevels[child] > dimToFound) {
                  std::vector<bool> dimFoundT(dimToFound, false);
                  childrenFinalT[n].clear();
                  childrenFinalT[n].resize(noNodes);
                  std::vector<bool> nodeDoneImplFind = nodeDoneImpl;
                  found = findStructGivenDim(child, dimToFound, false,
                                             nodeDoneImplFind, dimFoundT,
                                             childrenFinalT[n]);
                }
              } else
                found = true;
              if(found) {
                structOrigins[n] = child;
                nodeDoneImpl[child] = true;
                for(unsigned int j = 0; j < childrenFinalT[n].size(); ++j) {
                  for(auto &e : childrenFinalT[n][j]) {
                    childrenFinalOut[j].emplace_back(e);
                    nodeDoneImpl[e] = true;
                  }
                }
                break;
              }
            } // end for children[nodeOrigin]
            if(not found) {
              if(dimToFound <= 0) {
                structOrigins[n] = std::numeric_limits<ftm::idNode>::max();
                continue;
              }
              childrenFinalT[n].clear();
              childrenFinalT[n].resize(noNodes);
              createStructGivenDimRef(
                nodeOrigin, dimToFound, nodeDoneImpl, structOrigins[n],
                scalarsVectorOut, childrenFinalT[n], createStructGivenDimRef);
              for(unsigned int j = 0; j < childrenFinalT[n].size(); ++j) {
                for(auto &e : childrenFinalT[n][j]) {
                  if(e == structOrigins[n])
                    continue;
                  childrenFinalOut[j].emplace_back(e);
                }
              }
            }
          } // end for n
          // - Combine both structures
          if(structOrigins[0] == std::numeric_limits<ftm::idNode>::max()
             and structOrigins[1] == std::numeric_limits<ftm::idNode>::max()) {
            structOrigin = std::numeric_limits<ftm::idNode>::max();
            return;
          }
          bool firstIsParent = true;
          if(structOrigins[0] == std::numeric_limits<ftm::idNode>::max())
            firstIsParent = false;
          else if(structOrigins[1] == std::numeric_limits<ftm::idNode>::max())
            firstIsParent = true;
          else if(scalarsVectorOut[structOrigins[1] * 2 + 1]
                    - scalarsVectorOut[structOrigins[1] * 2]
                  > scalarsVectorOut[structOrigins[0] * 2 + 1]
                      - scalarsVectorOut[structOrigins[0] * 2])
            firstIsParent = false;
          structOrigin = (firstIsParent ? structOrigins[0] : structOrigins[1]);
          ftm::idNode modOrigin
            = (firstIsParent ? structOrigins[1] : structOrigins[0]);
          childrenFinalOut[nodeOrigin].emplace_back(structOrigin);
          if(modOrigin != std::numeric_limits<ftm::idNode>::max()) {
            childrenFinalOut[structOrigin].emplace_back(modOrigin);
            std::queue<std::array<ftm::idNode, 2>> queue;
            queue.emplace(std::array<ftm::idNode, 2>{modOrigin, structOrigin});
            while(!queue.empty()) {
              auto &nodeAndParent = queue.front();
              ftm::idNode node = nodeAndParent[0];
              ftm::idNode parent = nodeAndParent[1];
              queue.pop();
              adjustNestingScalars(scalarsVectorOut, node, parent);
              // Push children
              for(auto &child : childrenFinalOut[node])
                queue.emplace(std::array<ftm::idNode, 2>{child, node});
            }
          }
          return;
        };
      return createStructGivenDimImpl(
        _nodeOrigin, _dimToCreate, _nodeDone, _structOrigin, _scalarsVectorOut,
        _childrenFinalOut, createStructGivenDimImpl);
    };
  for(unsigned int i = 0; i < children.size(); ++i)
    sortChildren(i, nodeDone, children);
  Timer t_create;
  for(unsigned int i = 0; i < dimFound.size(); ++i) {
    if(dimFound[i])
      continue;
    ftm::idNode structOrigin;
    createStructGivenDim(
      startNode, i, nodeDone, structOrigin, scalarsVector, childrenFinal);
  }
}

void ttk::wae::printPairs(ftm::MergeTree<float> &mTree, bool useBD) {
  std::stringstream ss;
  if(mTree.tree.getRealNumberOfNodes() != 0)
    ss = mTree.tree.template printPairsFromTree<float>(useBD);
  else {
    std::vector<bool> nodeDone(mTree.tree.getNumberOfNodes(), false);
    for(unsigned int i = 0; i < mTree.tree.getNumberOfNodes(); ++i) {
      if(nodeDone[i])
        continue;
      std::tuple<ftm::idNode, ftm::idNode, float> pair
        = std::make_tuple(i, mTree.tree.getNode(i)->getOrigin(),
                          mTree.tree.getNodePersistence<float>(i));
      ss << std::get<0>(pair) << " ("
         << mTree.tree.getValue<float>(std::get<0>(pair)) << ") _ ";
      ss << std::get<1>(pair) << " ("
         << mTree.tree.getValue<float>(std::get<1>(pair)) << ") _ ";
      ss << std::get<2>(pair) << std::endl;
      nodeDone[i] = true;
      nodeDone[mTree.tree.getNode(i)->getOrigin()] = true;
    }
  }
  ss << std::endl;
  std::cout << ss.str();
}
