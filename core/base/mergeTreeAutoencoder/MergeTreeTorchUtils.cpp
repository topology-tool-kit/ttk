#include <MergeTreeTorchUtils.h>

using namespace std;
using namespace ttk;

#ifdef TTK_ENABLE_TORCH
using namespace torch::indexing;

void mtu::copyTensor(torch::Tensor &a, torch::Tensor &b) {
  b = a.detach().clone();
  b.requires_grad_(a.requires_grad());
}

void mtu::getDeltaProjTensor(torch::Tensor &diagTensor,
                             torch::Tensor &deltaProjTensor) {
  deltaProjTensor
    = (diagTensor.index({Slice(), 0}) + diagTensor.index({Slice(), 1})) / 2.0;
  deltaProjTensor = deltaProjTensor.reshape({-1, 1});
  deltaProjTensor = torch::cat({deltaProjTensor, deltaProjTensor}, 1);
}

void mtu::dataReorderingGivenMatching(mtu::TorchMergeTree<float> &tree,
                                      mtu::TorchMergeTree<float> &tree2,
                                      torch::Tensor &tree1ProjIndexer,
                                      torch::Tensor &tree2ReorderingIndexes,
                                      torch::Tensor &tree2ReorderedTensor,
                                      torch::Tensor &tree2DeltaProjTensor,
                                      torch::Tensor &tree1ReorderedTensor,
                                      torch::Tensor &tree2ProjIndexer,
                                      bool doubleReordering) {
  // Reorder tree2 tensor
  torch::Tensor tree2DiagTensor = tree2.tensor.reshape({-1, 2});
  tree2ReorderedTensor = torch::cat({tree2DiagTensor, torch::zeros({1, 2})});
  tree2ReorderedTensor = tree2ReorderedTensor.index({tree2ReorderingIndexes});

  // Create tree projection given matching
  torch::Tensor treeDiagTensor = tree.tensor.reshape({-1, 2});
  getDeltaProjTensor(treeDiagTensor, tree2DeltaProjTensor);
  tree2DeltaProjTensor = tree2DeltaProjTensor * tree1ProjIndexer;

  // Double reordering
  if(doubleReordering) {
    torch::Tensor tree1DeltaProjTensor;
    getDeltaProjTensor(tree2DiagTensor, tree1DeltaProjTensor);
    torch::Tensor tree2ProjIndexerR = tree2ProjIndexer.reshape({-1});
    tree1DeltaProjTensor = tree1DeltaProjTensor.index({tree2ProjIndexerR});
    tree1ReorderedTensor = torch::cat({treeDiagTensor, tree1DeltaProjTensor});
    tree1ReorderedTensor = tree1ReorderedTensor.reshape({-1, 1});
    torch::Tensor tree2UnmatchedTensor
      = tree2DiagTensor.index({tree2ProjIndexerR});
    tree2ReorderedTensor
      = torch::cat({tree2ReorderedTensor, tree2UnmatchedTensor});
    tree2DeltaProjTensor = torch::cat(
      {tree2DeltaProjTensor, torch::zeros_like(tree2UnmatchedTensor)});
  }

  // Reshape
  tree2ReorderedTensor = tree2ReorderedTensor.reshape({-1, 1});
  tree2DeltaProjTensor = tree2DeltaProjTensor.reshape({-1, 1});
}

void mtu::dataReorderingGivenMatching(mtu::TorchMergeTree<float> &tree,
                                      mtu::TorchMergeTree<float> &tree2,
                                      torch::Tensor &tree1ProjIndexer,
                                      torch::Tensor &tree2ReorderingIndexes,
                                      torch::Tensor &tree2ReorderedTensor,
                                      torch::Tensor &tree2DeltaProjTensor) {
  torch::Tensor tree1ReorderedTensor;
  torch::Tensor tree2ProjIndexer;
  bool doubleReordering = false;
  dataReorderingGivenMatching(tree, tree2, tree1ProjIndexer,
                              tree2ReorderingIndexes, tree2ReorderedTensor,
                              tree2DeltaProjTensor, tree1ReorderedTensor,
                              tree2ProjIndexer, doubleReordering);
}

void mtu::dataReorderingGivenMatching(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &tree2,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
  torch::Tensor &tree1ReorderedTensor,
  torch::Tensor &tree2ReorderedTensor,
  bool doubleReordering) {
  // Get tensor matching
  std::vector<int> tensorMatching;
  mtu::getTensorMatching(tree, tree2, matching, tensorMatching);
  torch::Tensor tree2ReorderingIndexes = torch::tensor(tensorMatching);
  torch::Tensor tree1ProjIndexer
    = (tree2ReorderingIndexes == -1).reshape({-1, 1});
  // Reorder tensor
  torch::Tensor tree2DeltaProjTensor;
  if(not doubleReordering) {
    dataReorderingGivenMatching(tree, tree2, tree1ProjIndexer,
                                tree2ReorderingIndexes, tree2ReorderedTensor,
                                tree2DeltaProjTensor);
  } else {
    std::vector<int> tensorMatching2;
    mtu::getInverseTensorMatching(tree, tree2, matching, tensorMatching);
    torch::Tensor tree1ReorderingIndexes = torch::tensor(tensorMatching);
    torch::Tensor tree2ProjIndexer
      = (tree1ReorderingIndexes == -1).reshape({-1, 1});
    dataReorderingGivenMatching(tree, tree2, tree1ProjIndexer,
                                tree2ReorderingIndexes, tree2ReorderedTensor,
                                tree2DeltaProjTensor, tree1ReorderedTensor,
                                tree2ProjIndexer, doubleReordering);
  }
  tree2ReorderedTensor = tree2ReorderedTensor + tree2DeltaProjTensor;
}

void mtu::dataReorderingGivenMatching(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &tree2,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
  torch::Tensor &tree2ReorderedTensor) {
  torch::Tensor tree1ReorderedTensor;
  bool doubleReordering = false;
  dataReorderingGivenMatching(tree, tree2, matching, tree1ReorderedTensor,
                              tree2ReorderedTensor, doubleReordering);
}

void mtu::meanBirthShift(torch::Tensor &diagTensor,
                         torch::Tensor &diagBaseTensor) {
  torch::Tensor birthShiftValue = diagBaseTensor.index({Slice(), 0}).mean()
                                  - diagTensor.index({Slice(), 0}).mean();
  torch::Tensor shiftTensor
    = torch::full({diagTensor.sizes()[0], 2}, birthShiftValue.item<float>());
  diagTensor.index_put_({None}, diagTensor + shiftTensor);
}

void mtu::meanBirthMaxPersShift(torch::Tensor &tensor,
                                torch::Tensor &baseTensor) {
  torch::Tensor diagTensor = tensor.reshape({-1, 2});
  torch::Tensor diagBaseTensor = baseTensor.reshape({-1, 2});
  // Shift to have same max pers
  torch::Tensor baseMaxPers
    = (diagBaseTensor.index({Slice(), 1}) - diagBaseTensor.index({Slice(), 0}))
        .max();
  torch::Tensor maxPers
    = (diagTensor.index({Slice(), 1}) - diagTensor.index({Slice(), 0})).max();
  torch::Tensor shiftTensor = (baseMaxPers - maxPers) / 2.0;
  shiftTensor = torch::stack({-shiftTensor, shiftTensor});
  diagTensor.index_put_({None}, diagTensor + shiftTensor);
  // Shift to have same birth mean
  meanBirthShift(diagTensor, diagBaseTensor);
}

void mtu::belowDiagonalPointsShift(torch::Tensor &tensor,
                                   torch::Tensor &backupTensor) {
  torch::Tensor oPDiag = tensor.reshape({-1, 2});
  torch::Tensor badPointsIndexer
    = (oPDiag.index({Slice(), 0}) > oPDiag.index({Slice(), 1}));
  torch::Tensor goodPoints = oPDiag.index({~badPointsIndexer});
  if(goodPoints.sizes()[0] == 0)
    goodPoints = backupTensor.reshape({-1, 2});
  torch::Tensor badPoints = oPDiag.index({badPointsIndexer});
  // Shift to be above diagonal with median pers
  torch::Tensor pers
    = (goodPoints.index({Slice(), 1}) - goodPoints.index({Slice(), 0}))
        .median();
  torch::Tensor shiftTensor
    = (torch::full({badPoints.sizes()[0], 1}, pers.item<float>())
       - badPoints.index({Slice(), 1}).reshape({-1, 1})
       + badPoints.index({Slice(), 0}).reshape({-1, 1}))
      / 2.0;
  shiftTensor = torch::cat({-shiftTensor, shiftTensor}, 1);
  badPoints = badPoints + shiftTensor;
  // Update tensor
  oPDiag.index_put_({badPointsIndexer}, badPoints);
  tensor = oPDiag.reshape({-1, 1}).detach();
}

void mtu::normalizeVectors(torch::Tensor &originTensor,
                           torch::Tensor &vectorsTensor) {
  torch::Tensor vSliced = vectorsTensor.index({Slice(2, None)});
  vSliced.index_put_({None}, vSliced / (originTensor[1] - originTensor[0]));
}

void mtu::normalizeVectors(mtu::TorchMergeTree<float> &origin,
                           std::vector<std::vector<double>> &vectors) {
  std::queue<ftm::idNode> queue;
  queue.emplace(origin.mTree.tree.getRoot());
  while(!queue.empty()) {
    ftm::idNode node = queue.front();
    queue.pop();
    if(not origin.mTree.tree.isRoot(node))
      for(unsigned int i = 0; i < 2; ++i)
        vectors[node][i] /= (origin.tensor[1] - origin.tensor[0]).item<float>();
    std::vector<ftm::idNode> children;
    origin.mTree.tree.getChildren(node, children);
    for(auto &child : children)
      queue.emplace(child);
  }
}

// Work only for persistence diagrams
bool mtu::isThereMissingPairs(mtu::TorchMergeTree<float> &interpolation) {
  float maxPers
    = interpolation.mTree.tree.template getMaximumPersistence<float>();
  torch::Tensor interTensor = interpolation.tensor;
  torch::Tensor indexer
    = torch::abs(interTensor.reshape({-1, 2}).index({Slice(), 0})
                 - interTensor.reshape({-1, 2}).index({Slice(), 1}))
      > (maxPers * 0.001 / 100.0);
  torch::Tensor indexed = interTensor.reshape({-1, 2}).index({indexer});
  return indexed.sizes()[0] > interpolation.mTree.tree.getRealNumberOfNodes();
}
#endif
