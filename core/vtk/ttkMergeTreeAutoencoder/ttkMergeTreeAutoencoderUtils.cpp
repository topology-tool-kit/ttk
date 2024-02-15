#include <MergeTreeAxesAlgorithmUtils.h>
#include <ttkMergeTreeAutoencoderUtils.h>

#ifdef TTK_ENABLE_TORCH
namespace ttk {
  namespace wae {
    void makeOneOutput(
      ttk::ftm::MergeTree<float> &tree,
      vtkUnstructuredGrid *treeNodes,
      std::vector<int> &treeNodeCorr,
      vtkDataSet *treeSegmentation,
      vtkSmartPointer<vtkUnstructuredGrid> &vtkOutputNode,
      vtkSmartPointer<vtkUnstructuredGrid> &vtkOutputArc,
      vtkSmartPointer<vtkDataSet> &vtkOutputSegmentation,
      unsigned int treeID,
      std::vector<std::tuple<std::string, std::vector<int>>> &customIntArrays,
      std::vector<std::tuple<std::string, std::vector<double>>>
        &customDoubleArrays,
      bool outputSegmentation,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      vtkOutputNode = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkOutputArc = vtkSmartPointer<vtkUnstructuredGrid>::New();

      ttkMergeTreeVisualization visuMakerBary;
      visuMakerBary.setShiftMode(-1); // Line
      visuMakerBary.setVtkOutputNode(vtkOutputNode);
      if(not isPersistenceDiagram)
        visuMakerBary.setVtkOutputArc(vtkOutputArc);
      else {
        visuMakerBary.setVtkOutputArc(vtkOutputNode);
        visuMakerBary.setIsPDSadMax(mixtureCoefficient == 0);
      }
      visuMakerBary.copyPointData(treeNodes, treeNodeCorr);
      for(auto &tup : customIntArrays)
        visuMakerBary.addCustomIntArray(std::get<0>(tup), std::get<1>(tup));
      for(auto &tup : customDoubleArrays)
        visuMakerBary.addCustomArray(std::get<0>(tup), std::get<1>(tup));
      visuMakerBary.setDebugLevel(debugLevel);
      visuMakerBary.setIsPersistenceDiagram(isPersistenceDiagram);
      visuMakerBary.setConvertedToDiagram(convertToDiagram);
      if(treeNodes) {
        visuMakerBary.setTreesNodes(treeNodes);
        visuMakerBary.setTreesNodeCorrMesh(treeNodeCorr);
      }
      if(outputSegmentation) {
        vtkOutputSegmentation = vtkSmartPointer<vtkUnstructuredGrid>::New();
        visuMakerBary.setTreesSegmentation(treeSegmentation);
        visuMakerBary.setPlanarLayout(false);
        visuMakerBary.setOutputSegmentation(true);
        visuMakerBary.setVtkOutputSegmentation(vtkOutputSegmentation);
      } else {
        visuMakerBary.setPlanarLayout(true);
      }
      visuMakerBary.setISampleOffset(treeID);
      visuMakerBary.setOutputTreeNodeId(true);
      visuMakerBary.makeTreesOutput<float>(&(tree.tree));
    }

    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      std::vector<vtkUnstructuredGrid *> &treesNodesT,
      std::vector<std::vector<int>> &treesNodeCorr,
      std::vector<vtkDataSet *> &treesSegmentationT,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      vtkSmartPointer<vtkMultiBlockDataSet> allNodes
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkSmartPointer<vtkMultiBlockDataSet> allArcs;
      if(not isPersistenceDiagram) {
        allArcs = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      }
      bool outputSegmentation
        = !treesSegmentationT.empty() and treesSegmentationT[0];
      vtkSmartPointer<vtkMultiBlockDataSet> allSegs;
      if(outputSegmentation) {
        allSegs = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      }
      int shift = 0;
      for(unsigned int i = 0; i < trees.size(); ++i) {
        if(trees[i]->tree.template getMaximumPersistence<float>() == 0) {
          ++shift;
          continue;
        }
        vtkUnstructuredGrid *treeNodes = nullptr;
        vtkDataSet *treeSegmentation = nullptr;
        std::vector<int> treeNodeCorr;
        if(outputSegmentation) {
          treeSegmentation = treesSegmentationT[i];
        }
        if(i < treesNodesT.size()) {
          treeNodes = treesNodesT[i];
          treeNodeCorr = treesNodeCorr[i];
        }
        vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode, vtkOutputArc;
        vtkSmartPointer<vtkDataSet> vtkOutputSegmentation;
        makeOneOutput(*(trees[i]), treeNodes, treeNodeCorr, treeSegmentation,
                      vtkOutputNode, vtkOutputArc, vtkOutputSegmentation, i,
                      customIntArrays[i], customDoubleArrays[i],
                      outputSegmentation, mixtureCoefficient,
                      isPersistenceDiagram, convertToDiagram, debugLevel);
        allNodes->SetBlock(i - shift, vtkOutputNode);
        if(not isPersistenceDiagram)
          allArcs->SetBlock(i - shift, vtkOutputArc);
        if(outputSegmentation)
          allSegs->SetBlock(i - shift, vtkOutputSegmentation);
      }
      if(not isPersistenceDiagram) {
        output->SetNumberOfBlocks(2);
        output->SetBlock(0, allNodes);
        output->SetBlock(1, allArcs);
        if(outputSegmentation)
          output->SetBlock(2, allSegs);
      } else {
        if(not outputSegmentation) {
          output->ShallowCopy(allNodes);
        } else {
          output->SetNumberOfBlocks(2);
          output->SetBlock(0, allNodes);
          output->SetBlock(1, allSegs);
        }
      }
    }

    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      std::vector<vtkUnstructuredGrid *> &treesNodesT,
      std::vector<std::vector<int>> &treesNodeCorr,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      std::vector<vtkDataSet *> treesSegmentationT;
      makeManyOutput(trees, treesNodesT, treesNodeCorr, treesSegmentationT,
                     output, customIntArrays, customDoubleArrays,
                     mixtureCoefficient, isPersistenceDiagram, convertToDiagram,
                     debugLevel);
    }

    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      std::vector<vtkUnstructuredGrid *> treesNodesT;
      std::vector<std::vector<int>> treesNodeCorr;
      makeManyOutput(trees, treesNodesT, treesNodeCorr, output, customIntArrays,
                     customDoubleArrays, mixtureCoefficient,
                     isPersistenceDiagram, convertToDiagram, debugLevel);
    }

    void makeManyOutput(std::vector<ttk::ftm::MergeTree<float> *> &trees,
                        vtkSmartPointer<vtkMultiBlockDataSet> &output,
                        double mixtureCoefficient,
                        bool isPersistenceDiagram,
                        bool convertToDiagram,
                        int debugLevel) {
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        customIntArrays(trees.size());
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        customDoubleArrays(trees.size());
      makeManyOutput(trees, output, customIntArrays, customDoubleArrays,
                     mixtureCoefficient, isPersistenceDiagram, convertToDiagram,
                     debugLevel);
    }

    void computeTrackingInformation(
      std::vector<mtu::TorchMergeTree<float>> &origins,
      std::vector<mtu::TorchMergeTree<float>> &originsPrime,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invOriginsMatchingVectorT,
      bool isPersistenceDiagram,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<double> &originPersPercent,
      std::vector<double> &originPersDiff,
      std::vector<int> &originPersistenceOrder) {
      unsigned int originsMatchingSize = originsMatchingVectorT.size();
      originsMatchingVector.resize(originsMatchingSize);
      originsPersPercent.resize(originsMatchingSize);
      originsPersDiff.resize(originsMatchingSize);
      for(unsigned int l = 0; l < originsMatchingSize; ++l) {
        auto &tree2 = (l == 0 ? originsPrime[0] : originsPrime[l]);
        originsMatchingVector[l] = invOriginsMatchingVectorT[l];
        if(l != 0) {
          for(unsigned int i = 0; i < originsMatchingVector[l].size(); ++i)
            if(originsMatchingVector[l][i]
               < originsMatchingVector[l - 1].size())
              originsMatchingVector[l][i]
                = originsMatchingVector[l - 1][originsMatchingVector[l][i]];
        }
        originsPersPercent[l].resize(tree2.mTree.tree.getNumberOfNodes());
        originsPersDiff[l].resize(tree2.mTree.tree.getNumberOfNodes());
        for(unsigned int i = 0; i < originsMatchingVector[l].size(); ++i) {
          if(originsMatchingVector[l][i]
             >= origins[0].mTree.tree.getNumberOfNodes())
            continue;
          auto pers = origins[0].mTree.tree.template getNodePersistence<float>(
            originsMatchingVector[l][i]);
          auto treePers
            = tree2.mTree.tree.template getNodePersistence<float>(i);
          originsPersPercent[l][i] = treePers * 100 / pers;
          originsPersDiff[l][i] = treePers - pers;
        }
      }

      originPersPercent.resize(origins[0].mTree.tree.getNumberOfNodes());
      originPersDiff.resize(origins[0].mTree.tree.getNumberOfNodes());
      std::vector<ttk::ftm::idNode> originMatchingVector;
      for(unsigned int l = 0; l < originsMatchingSize; ++l) {
        std::vector<ttk::ftm::idNode> &originMatchingVectorT
          = originsMatchingVectorT[l];
        if(l == 0) {
          originMatchingVector = originMatchingVectorT;
        } else {
          for(unsigned int i = 0; i < originMatchingVector.size(); ++i)
            if(originMatchingVector[i] < originMatchingVectorT.size())
              originMatchingVector[i]
                = originMatchingVectorT[originMatchingVector[i]];
        }
      }
      unsigned int l2 = originsMatchingSize - 1;
      for(unsigned int i = 0; i < originMatchingVector.size(); ++i) {
        if(originMatchingVector[i] < originsPersDiff[l2].size()) {
          originPersPercent[i]
            = originsPersPercent[l2][originMatchingVector[i]];
          originPersDiff[i] = originsPersDiff[l2][originMatchingVector[i]];
        }
      }

      originPersistenceOrder.resize(
        origins[0].mTree.tree.getNumberOfNodes(), -1);
      std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, float>>
        pairsBary;
      bool useBD = isPersistenceDiagram;
      origins[0].mTree.tree.template getPersistencePairsFromTree<float>(
        pairsBary, useBD);
      for(unsigned int j = 0; j < pairsBary.size(); ++j) {
        int index = pairsBary.size() - 1 - j;
        originPersistenceOrder[std::get<0>(pairsBary[j])] = index;
        originPersistenceOrder[std::get<1>(pairsBary[j])] = index;
      }
    }

    void computeCustomArrays(
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
      std::vector<std::vector<double>> &persCorrelationMatrix,
      std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
        &invDataMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invReconstMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<int> &originPersistenceOrder,
      unsigned int l,
      unsigned int lShift,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays) {
      unsigned int originsMatchingSize = originsMatchingVectorT.size();
      unsigned int dataMatchingSize = invDataMatchingVectorT.size();
      unsigned int lShifted = l + lShift;
      std::vector<std::vector<ttk::ftm::idNode>> matchingVectors(recs.size());
      std::vector<std::vector<double>> dataPersPercent, dataPersDiff;
      std::vector<std::vector<int>> dataOriginPersOrder;
      std::vector<std::vector<std::vector<double>>> dataCorrelation;
      if(lShifted < dataMatchingSize) {
        for(unsigned int i = 0; i < recs.size(); ++i) {
          matchingVectors[i] = invDataMatchingVectorT[lShifted][i];
          if(lShifted != 0) {
            for(unsigned int j = 0; j < matchingVectors[i].size(); ++j)
              if(matchingVectors[i][j]
                 < originsMatchingVector[lShifted - 1].size())
                matchingVectors[i][j]
                  = originsMatchingVector[lShifted - 1][matchingVectors[i][j]];
          }
        }
      }
      if(lShifted == 0 or lShifted == dataMatchingSize - 1) {
        dataPersPercent.resize(recs.size());
        dataPersDiff.resize(recs.size());
        for(unsigned int i = 0; i < recs.size(); ++i) {
          dataPersPercent[i].resize(recs[i][l].mTree.tree.getNumberOfNodes());
          dataPersDiff[i].resize(recs[i][l].mTree.tree.getNumberOfNodes());
          std::vector<ttk::ftm::idNode> matchingVector;
          if(lShifted == 0) {
            matchingVector = matchingVectors[i];
            for(unsigned int l2 = 0; l2 < originsMatchingSize; ++l2) {
              std::vector<ttk::ftm::idNode> &originMatchingVector
                = originsMatchingVectorT[l2];
              for(unsigned int j = 0; j < matchingVector.size(); ++j)
                if(matchingVector[j] < originMatchingVector.size())
                  matchingVector[j] = originMatchingVector[matchingVector[j]];
            }
          } else {
            matchingVector = invDataMatchingVectorT[lShifted][i];
          }
          unsigned int l2 = originsMatchingSize - 1;
          for(unsigned int j = 0; j < matchingVector.size(); ++j) {
            if(matchingVector[j] < originsPersDiff[l2].size()) {
              dataPersDiff[i][j] = originsPersDiff[l2][matchingVector[j]];
              dataPersPercent[i][j] = originsPersPercent[l2][matchingVector[j]];
            }
          }
        }

        if(lShifted == 0) {
          dataCorrelation.resize(recs.size());
          for(unsigned int i = 0; i < recs.size(); ++i) {
            dataCorrelation[i].resize(persCorrelationMatrix[0].size());
            for(unsigned int j = 0; j < persCorrelationMatrix[0].size(); ++j) {
              dataCorrelation[i][j].resize(
                recs[i][l].mTree.tree.getNumberOfNodes());
              for(unsigned int k = 0; k < matchingVectors[i].size(); ++k) {
                if(matchingVectors[i][k] < persCorrelationMatrix.size())
                  dataCorrelation[i][j][k]
                    = persCorrelationMatrix[matchingVectors[i][k]][j];
              }
            }
          }
        }
      }

      if(lShifted == 0 or lShifted == dataMatchingSize - 1
         or l == recs[0].size() - 1) {
        dataOriginPersOrder.resize(recs.size());
        for(unsigned int i = 0; i < recs.size(); ++i) {
          std::vector<ttk::ftm::idNode> &matchingVector = matchingVectors[i];
          if(l == recs[0].size() - 1) {
            matchingVector = invReconstMatchingVectorT[i];
            std::vector<ttk::ftm::idNode> &matchingVectorT
              = invDataMatchingVectorT[0][i];
            for(unsigned int j = 0; j < matchingVector.size(); ++j)
              if(matchingVector[j] < matchingVectorT.size())
                matchingVector[j] = matchingVectorT[matchingVector[j]];
          }
          dataOriginPersOrder[i].resize(
            recs[i][l].mTree.tree.getNumberOfNodes());
          for(unsigned int j = 0; j < matchingVector.size(); ++j) {
            if(matchingVector[j] < originPersistenceOrder.size())
              dataOriginPersOrder[i][j]
                = originPersistenceOrder[matchingVector[j]];
            else
              dataOriginPersOrder[i][j] = -1;
          }
        }
      }

      for(unsigned int i = 0; i < recs.size(); ++i) {
        if(lShifted < dataMatchingSize) {
          std::vector<int> customArrayMatching;
          for(auto &e : matchingVectors[i])
            customArrayMatching.emplace_back(e);
          std::string name{"OriginTrueNodeId"};
          customIntArrays[i].emplace_back(
            std::make_tuple(name, customArrayMatching));
          if(lShifted == 0 or lShifted == dataMatchingSize - 1) {
            std::string name2{"OriginPersPercent"};
            customDoubleArrays[i].emplace_back(
              std::make_tuple(name2, dataPersPercent[i]));
            std::string name3{"OriginPersDiff"};
            customDoubleArrays[i].emplace_back(
              std::make_tuple(name3, dataPersDiff[i]));
          }
          if(lShifted == 0) {
            for(unsigned int j = 0; j < dataCorrelation[i].size(); ++j) {
              std::string name2 = ttk::axa::getTableCorrelationPersName(
                dataCorrelation[i].size(), j);
              customDoubleArrays[i].emplace_back(
                std::make_tuple(name2, dataCorrelation[i][j]));
            }
          }
        }
        if(lShifted == 0 or lShifted == dataMatchingSize - 1
           or l == recs[0].size() - 1) {
          std::string name4{"OriginPersOrder"};
          customIntArrays[i].emplace_back(
            std::make_tuple(name4, dataOriginPersOrder[i]));
        }
      }
    }
  } // namespace wae
} // namespace ttk
#endif
