/// \ingroup base
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// \brief Provide methods related to Merge Trees and their management with
/// Torch.

#pragma once

#include <Debug.h>
#include <FTMTreeUtils.h>
#include <FTMTree_MT.h>
#include <Geometry.h>
#include <MergeTreeUtils.h>

#ifdef TTK_ENABLE_TORCH
#include <torch/torch.h>
#endif

namespace ttk {

  namespace mtu {

#ifdef TTK_ENABLE_TORCH
    /**
     * @brief Copy tensor.
     *
     * @param[in] a input tensor.
     * @param[out] b copied output tensor.
     */
    void copyTensor(torch::Tensor &a, torch::Tensor &b);

    template <typename dataType>
    struct TorchMergeTree {
      ftm::MergeTree<dataType> mTree;
      torch::Tensor tensor;
      std::vector<unsigned int> nodeCorr;
      std::vector<ftm::idNode> parentsOri;
    };

    /**
     * @brief Get a tensor of the pairs projected to the diagonal.
     *
     * @param[in] diagTensor tensor of persistence pairs.
     * @param[out] deltaProjTensor tensor of projected pairs.
     */
    void getDeltaProjTensor(torch::Tensor &diagTensor,
                            torch::Tensor &deltaProjTensor);

    /**
     * @brief Compute the reordering of torch tensors given a matching.
     *
     * @param[in] tree first torch merge tree.
     * @param[in] tree2 second torch merge tree.
     * @param[in] tree1ProjIndexer torch tensor indexing the destroyed pairs in
     * first tree.
     * @param[in] tree2ReorderingIndexes tensor reordering the second tree given
     * a matching.
     * @param[out] tree2ReorderedTensor reordered torch tensor of second tree
     * (with zero on non-matched pairs).
     * @param[out] tree2DeltaProjTensor tensor of the projected pairs on the
     * diagonal of the first tree in the second tree.
     * @param[out] tree1ReorderedTensor reordered torch tensor of first tree.
     * @param[in] tree2ProjIndexer torch tensor indexing the destroyed pairs in
     * second tree.
     * @param[in] doubleReordering choose to also reorder first tree.
     */
    void dataReorderingGivenMatching(mtu::TorchMergeTree<float> &tree,
                                     mtu::TorchMergeTree<float> &tree2,
                                     torch::Tensor &tree1ProjIndexer,
                                     torch::Tensor &tree2ReorderingIndexes,
                                     torch::Tensor &tree2ReorderedTensor,
                                     torch::Tensor &tree2DeltaProjTensor,
                                     torch::Tensor &tree1ReorderedTensor,
                                     torch::Tensor &tree2ProjIndexer,
                                     bool doubleReordering = true);

    /**
     * @brief Compute the reordering of torch tensors given a matching.
     *
     * @param[in] tree first torch merge tree.
     * @param[in] tree2 second torch merge tree.
     * @param[in] tree1ProjIndexer torch tensor indexing the destroyed pairs in
     * first tree.
     * @param[in] tree2ReorderingIndexes tensor reordering the second tree given
     * a matching.
     * @param[out] tree2ReorderedTensor reordered torch tensor of second tree
     * (with zero on non-matched pairs).
     * @param[out] tree2DeltaProjTensor tensor of the projected pairs on the
     * diagonal of the first tree in the second tree.
     */
    void dataReorderingGivenMatching(mtu::TorchMergeTree<float> &tree,
                                     mtu::TorchMergeTree<float> &tree2,
                                     torch::Tensor &tree1ProjIndexer,
                                     torch::Tensor &tree2ReorderingIndexes,
                                     torch::Tensor &tree2ReorderedTensor,
                                     torch::Tensor &tree2DeltaProjTensor);

    /**
     * @brief Compute the reordering of torch tensors given a matching.
     *
     * @param[in] tree first torch merge tree.
     * @param[in] tree2 second torch merge tree.
     * @param[in] matching matching between trees.
     * @param[out] tree1ReorderedTensor reordered torch tensor of first tree.
     * @param[out] tree2ReorderedTensor reordered torch tensor of second tree.
     * @param[in] doubleReordering choose to also reorder first tree.
     */
    void dataReorderingGivenMatching(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      torch::Tensor &tree1ReorderedTensor,
      torch::Tensor &tree2ReorderedTensor,
      bool doubleReordering = true);

    /**
     * @brief Compute the reordering of a torch tensor given a matching.
     *
     * @param[in] tree first torch merge tree.
     * @param[in] tree2 second torch merge tree.
     * @param[in] matching matching between trees.
     * @param[out] tree2ReorderedTensor reordered torch tensor of second tree.
     */
    void dataReorderingGivenMatching(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      torch::Tensor &tree2ReorderedTensor);

    /**
     * @brief Shift pairs to have same birth mean than a reference tensor.
     *
     * @param[in] diagTensor tensor of persistence pairs.
     * @param[out] diagBaseTensor reference tensor.
     */
    void meanBirthShift(torch::Tensor &diagTensor,
                        torch::Tensor &diagBaseTensor);

    /**
     * @brief Shift pairs to have same max persistence and same birth mean than
     * a reference tensor.
     *
     * @param[in] tensor tensor of persistence pairs.
     * @param[out] baseTensor reference tensor.
     */
    void meanBirthMaxPersShift(torch::Tensor &tensor,
                               torch::Tensor &baseTensor);

    /**
     * @brief Shift pairs below diagonal above it with a median persistence.
     *
     * @param[in] tensor tensor of persistence pairs.
     * @param[out] backupTensor reference tensor if all pairs are below
     * diagonal.
     */
    void belowDiagonalPointsShift(torch::Tensor &tensor,
                                  torch::Tensor &backupTensor);

    /**
     * @brief Normalize an axis vector.
     *
     * @param[in] originTensor input torch tensor.
     * @param[in] vectorsTensor vectors tensor to be normalized.
     */
    void normalizeVectors(torch::Tensor &originTensor,
                          torch::Tensor &vectorsTensor);

    /**
     * @brief Normalize axis vectors.
     *
     * @param[in] origin input torch merge tree.
     * @param[in] vectors vectors to be normalized.
     */
    void normalizeVectors(mtu::TorchMergeTree<float> &origin,
                          std::vector<std::vector<double>> &vectors);

    /**
     * @brief Get the latent layer index.
     *
     * @return the latent layer index.
     */
    unsigned int getLatentLayerIndex();

    /**
     * @brief Test if pairs are missing (testing function).
     *
     * @return return true if pairs are missing.
     */
    bool isThereMissingPairs(mtu::TorchMergeTree<float> &interpolation);

    /**
     * @brief Copy a torch merge tree.
     *
     * @param[in] tmTree input torch merge tree.
     * @param[out] out output copied torch merge tree.
     */
    template <class dataType>
    void copyTorchMergeTree(TorchMergeTree<dataType> &tmTree,
                            TorchMergeTree<dataType> &out) {
      out.mTree = ftm::copyMergeTree<dataType>(tmTree.mTree);
      copyTensor(tmTree.tensor, out.tensor);
      out.nodeCorr = tmTree.nodeCorr;
      out.parentsOri = tmTree.parentsOri;
    }

    /**
     * @brief Convert a merge tree to a torch tensor.
     *
     * @param[in] mTree input merge tree.
     * @param[out] tensor output tensor.
     * @param[out] nodeCorr correspondence between the nodes of the tree and the
     * row index in the tensor.
     * @param[in] normalize if normalization is used.
     * @param[in] revNodeCorr pointer to enforce a specific correspondence
     * between nodes and row indices.
     * @param[in] revNodeCorrSize size of the revNodeCorr array.
     */
    template <class dataType>
    void mergeTreeToTorchTensor(ftm::MergeTree<dataType> &mTree,
                                torch::Tensor &tensor,
                                std::vector<unsigned int> &nodeCorr,
                                bool normalize,
                                unsigned int *revNodeCorr = nullptr,
                                unsigned int revNodeCorrSize = 0) {
      nodeCorr.clear();
      nodeCorr.assign(mTree.tree.getNumberOfNodes(),
                      std::numeric_limits<unsigned int>::max());
      std::vector<dataType> scalars;
      std::queue<ftm::idNode> queue;
      if(revNodeCorrSize == 0)
        queue.emplace(mTree.tree.getRoot());
      else {
        for(unsigned int i = 0; i < revNodeCorrSize; i += 2)
          queue.emplace(revNodeCorr[i]);
      }
      unsigned int cpt = 0;
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        auto birthDeath = ttk::getParametrizedBirthDeath<dataType>(
          &(mTree.tree), node, normalize);
        auto birth = std::get<0>(birthDeath);
        auto death = std::get<1>(birthDeath);
        scalars.emplace_back(birth);
        scalars.emplace_back(death);

        nodeCorr[node] = cpt;
        ++cpt;

        if(revNodeCorrSize == 0) {
          std::vector<ftm::idNode> children;
          mTree.tree.getChildren(node, children);
          for(auto child : children)
            queue.emplace(child);
        }
      }
      tensor = torch::tensor(scalars).reshape({-1, 1});
    }

    /**
     * @brief Convert a merge tree to a torch tensor.
     *
     * @param[in] mTree input merge tree.
     * @param[out] tensor output tensor.
     * @param[in] normalize if normalization is used.
     */
    template <class dataType>
    void mergeTreeToTorchTensor(ftm::MergeTree<dataType> &mTree,
                                torch::Tensor &tensor,
                                bool normalize) {
      std::vector<unsigned int> nodeCorr;
      mergeTreeToTorchTensor<dataType>(mTree, tensor, nodeCorr, normalize);
    }

    /**
     * @brief Get a vector assigning to each node the index of its parent.
     *
     * @param[in] mTree input merge tree.
     * @param[out] parents parents vector.
     */
    template <class dataType>
    void getParentsVector(ftm::MergeTree<dataType> &mTree,
                          std::vector<ftm::idNode> &parents) {
      parents.resize(mTree.tree.getNumberOfNodes());
      std::fill(parents.begin(), parents.end(),
                std::numeric_limits<ftm::idNode>::max());
      std::queue<ftm::idNode> queue;
      queue.emplace(mTree.tree.getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        if(!mTree.tree.isRoot(node))
          parents[node] = mTree.tree.getParentSafe(node);
        std::vector<ftm::idNode> children;
        mTree.tree.getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
    }

    /**
     * @brief Convert a merge tree to a torch merge tree.
     *
     * @param[in] mTree input merge tree.
     * @param[out] out output torch merge tree.
     * @param[in] normalize if normalization is used.
     */
    template <class dataType>
    void mergeTreeToTorchTree(ftm::MergeTree<dataType> &mTree,
                              TorchMergeTree<dataType> &out,
                              bool normalize) {
      out.mTree = copyMergeTree(mTree);
      getParentsVector(out.mTree, out.parentsOri);
      mergeTreeToTorchTensor<dataType>(
        out.mTree, out.tensor, out.nodeCorr, normalize);
    }

    /**
     * @brief Convert merge trees to torch merge trees.
     *
     * @param[in] mTrees input merge trees.
     * @param[out] out output torch merge trees.
     * @param[in] normalize if normalization is used.
     */
    template <class dataType>
    void mergeTreesToTorchTrees(std::vector<ftm::MergeTree<dataType>> &mTrees,
                                std::vector<TorchMergeTree<dataType>> &out,
                                bool normalize) {
      out.resize(mTrees.size());
      for(unsigned int i = 0; i < mTrees.size(); ++i)
        mergeTreeToTorchTree<dataType>(mTrees[i], out[i], normalize);
    }

    /**
     * @brief Convert a merge tree to a torch merge tree.
     *
     * @param[in] mTree input merge tree.
     * @param[out] out output torch merge tree.
     * @param[in] normalize if normalization is used.
     * @param[in] revNodeCorr pointer to enforce a specific correspondence
     * between nodes and row indices.
     * @param[in] revNodeCorrSize size of the revNodeCorr array.
     */
    template <class dataType>
    void mergeTreeToTorchTree(ftm::MergeTree<dataType> &mTree,
                              TorchMergeTree<dataType> &out,
                              bool normalize,
                              unsigned int *revNodeCorr,
                              unsigned int revNodeCorrSize) {
      out.mTree = copyMergeTree(mTree);
      getParentsVector(out.mTree, out.parentsOri);
      mergeTreeToTorchTensor<dataType>(out.mTree, out.tensor, out.nodeCorr,
                                       normalize, revNodeCorr, revNodeCorrSize);
    }

    /**
     * @brief Convert merge trees to torch merge trees.
     *
     * @param[in] mTrees input merge trees.
     * @param[out] out output torch merge trees.
     * @param[in] normalize if normalization is used.
     * @param[in] revNodeCorr pointer to enforce a specific correspondence
     * between nodes and row indices.
     * @param[in] revNodeCorrSize size of the revNodeCorr array.
     */
    template <class dataType>
    void mergeTreesToTorchTrees(std::vector<ftm::MergeTree<dataType>> &mTrees,
                                std::vector<TorchMergeTree<dataType>> &out,
                                bool normalize,
                                std::vector<unsigned int *> &allRevNodeCorr,
                                std::vector<unsigned int> &allRevNodeCorrSize) {
      out.resize(mTrees.size());
      for(unsigned int i = 0; i < mTrees.size(); ++i)
        mergeTreeToTorchTree<dataType>(mTrees[i], out[i], normalize,
                                       allRevNodeCorr[i],
                                       allRevNodeCorrSize[i]);
    }

    /**
     * @brief Convert a torch tensor to a merge tree.
     *
     * @param[in] tmt input torch merge tree.
     * @param[in] normalized if normalization is used.
     * @param[out] mTreeOut output merge tree
     */
    template <class dataType>
    bool torchTensorToMergeTree(TorchMergeTree<dataType> &tmt,
                                bool normalized,
                                ftm::MergeTree<dataType> &mTreeOut) {
      std::vector<unsigned int> &nodeCorr = tmt.nodeCorr;
      torch::Tensor &tensor = tmt.tensor;
      std::vector<ftm::idNode> &parentsOri = tmt.parentsOri;

      mTreeOut = ttk::ftm::copyMergeTree<dataType>(tmt.mTree);
      if(parentsOri.empty())
        std::cout << "[torchTensorToMergeTree] parentsOri.empty()" << std::endl;
      for(unsigned int i = 0; i < parentsOri.size(); ++i)
        if(parentsOri[i] < mTreeOut.tree.getNumberOfNodes()
           and mTreeOut.tree.isRoot(i))
          mTreeOut.tree.setParent(i, parentsOri[i]);
      ftm::idNode root = mTreeOut.tree.getRoot();
      if(root >= mTreeOut.tree.getNumberOfNodes())
        return true;

      bool isJT = tmt.mTree.tree.template isJoinTree<dataType>();
      std::vector<dataType> tensorVec(
        tensor.data_ptr<float>(), tensor.data_ptr<float>() + tensor.numel());
      std::vector<dataType> scalarsVector;
      ttk::ftm::getTreeScalars<dataType>(mTreeOut, scalarsVector);
      std::queue<ftm::idNode> queue;
      queue.emplace(root);
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        auto birthValue = tensorVec[nodeCorr[node] * 2];
        auto deathValue = tensorVec[nodeCorr[node] * 2 + 1];
        if(normalized and !mTreeOut.tree.isRoot(node)) {
          ftm::idNode nodeParent = mTreeOut.tree.getParentSafe(node);
          ftm::idNode nodeParentOrigin
            = mTreeOut.tree.getNode(nodeParent)->getOrigin();
          ftm::idNode birthParentNode = (isJT ? nodeParentOrigin : nodeParent);
          ftm::idNode deathParentNode = (isJT ? nodeParent : nodeParentOrigin);
          auto birthParent = scalarsVector[birthParentNode];
          auto deathParent = scalarsVector[deathParentNode];
          birthValue = birthValue * (deathParent - birthParent) + birthParent;
          deathValue = deathValue * (deathParent - birthParent) + birthParent;
        }
        ftm::idNode nodeOrigin = mTreeOut.tree.getNode(node)->getOrigin();
        scalarsVector[node] = (isJT ? deathValue : birthValue);
        scalarsVector[nodeOrigin] = (isJT ? birthValue : deathValue);

        std::vector<ftm::idNode> children;
        mTreeOut.tree.getChildren(node, children);
        for(auto &child : children)
          queue.emplace(child);
      }
      ftm::setTreeScalars<dataType>(mTreeOut, scalarsVector);
      return false;
    }

    /**
     * @brief Assign the parent of each node of a torch merge tree given its
     * parents vecotr.
     *
     * @param[in] tmt input torch merge tree.
     */
    template <class dataType>
    void fillMergeTreeStructure(TorchMergeTree<dataType> &tmt) {
      std::vector<ftm::idNode> &parentsOri = tmt.parentsOri;
      for(unsigned int i = 0; i < parentsOri.size(); ++i)
        if(parentsOri[i] < tmt.mTree.tree.getNumberOfNodes()
           and tmt.mTree.tree.isRoot(i))
          tmt.mTree.tree.setParent(i, parentsOri[i]);
    }

    /**
     * @brief Create the reverse correspondence between the row indices of the
     * tensor a torch merge tree to its nodes.
     *
     * @param[in] tmt input torch merge tree.
     * @param[out] revNodeCorr output node correspondence.
     */
    template <class dataType>
    void getReverseTorchNodeCorr(TorchMergeTree<dataType> &tmt,
                                 std::vector<unsigned int> &revNodeCorr) {
      revNodeCorr.clear();
      revNodeCorr.assign(
        tmt.tensor.sizes()[0], std::numeric_limits<unsigned int>::max());
      for(unsigned int i = 0; i < tmt.nodeCorr.size(); ++i) {
        if(tmt.nodeCorr[i] != std::numeric_limits<unsigned int>::max()) {
          revNodeCorr[tmt.nodeCorr[i] * 2] = i;
          revNodeCorr[tmt.nodeCorr[i] * 2 + 1] = i;
        }
      }
    }

    /**
     * @brief Convert an axis vector to a torch tensor.
     *
     * @param[in] mTree input merge tree for the tensor ordering.
     * @param[in] v axis vector.
     * @param[out] tensor output tensor.
     */
    template <class dataType>
    void axisVectorToTorchTensor(ftm::MergeTree<dataType> &mTree,
                                 std::vector<std::vector<double>> &v,
                                 torch::Tensor &tensor) {
      std::vector<double> v_flatten;
      std::queue<ftm::idNode> queue;
      queue.emplace(mTree.tree.getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        v_flatten.emplace_back(v[node][0]);
        v_flatten.emplace_back(v[node][1]);

        std::vector<ftm::idNode> children;
        mTree.tree.getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
      tensor = torch::tensor(v_flatten);
      tensor = tensor.reshape({-1, 1});
    }

    /**
     * @brief Convert axis vectors to torch tensors.
     *
     * @param[in] mTree input merge tree for the tensor ordering.
     * @param[in] vS axis vectors.
     * @param[out] tensor output tensors.
     */
    template <class dataType>
    void axisVectorsToTorchTensor(
      ftm::MergeTree<dataType> &mTree,
      std::vector<std::vector<std::vector<double>>> &vS,
      torch::Tensor &tensor) {
      std::vector<torch::Tensor> allTensors;
      for(auto &v : vS) {
        torch::Tensor t;
        axisVectorToTorchTensor(mTree, v, t);
        allTensors.emplace_back(t);
      }
      tensor = torch::cat(allTensors, 1);
    }

    /**
     * @brief Create a matching tensor.
     *
     * @param[in] a first input torch merge tree.
     * @param[in] b second input torch merge tree.
     * @param[in] matching matching between input merge trees.
     * @param[out] tensorMatching ouput tensor matching.
     */
    template <class dataType>
    void getTensorMatching(
      TorchMergeTree<dataType> &a,
      TorchMergeTree<dataType> &b,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<int> &tensorMatching) {
      tensorMatching.clear();
      tensorMatching.resize((int)(a.tensor.sizes()[0] / 2), -1);
      for(auto &match : matching) {
        auto match1 = std::get<0>(match);
        auto match2 = std::get<1>(match);
        if(a.nodeCorr[match1] < tensorMatching.size())
          tensorMatching[a.nodeCorr[match1]] = b.nodeCorr[match2];
      }
    }

    /**
     * @brief Create an inverse matching tensor.
     *
     * @param[in] a first input torch merge tree.
     * @param[in] b second input torch merge tree.
     * @param[in] matching matching between input merge trees.
     * @param[out] tensorMatching ouput tensor matching.
     */
    template <class dataType>
    void getInverseTensorMatching(
      TorchMergeTree<dataType> &a,
      TorchMergeTree<dataType> &b,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<int> &tensorMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> invMatching(
        matching.size());
      for(unsigned int i = 0; i < matching.size(); ++i)
        invMatching[i]
          = std::make_tuple(std::get<1>(matching[i]), std::get<0>(matching[i]),
                            std::get<2>(matching[i]));
      getTensorMatching(b, a, invMatching, tensorMatching);
    }
#endif

  } // namespace mtu

} // namespace ttk
