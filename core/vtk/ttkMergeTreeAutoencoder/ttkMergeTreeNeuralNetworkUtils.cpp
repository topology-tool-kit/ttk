#include <MergeTreeAxesAlgorithmUtils.h>
#include <ttkMergeTreeNeuralNetworkUtils.h>

#include <vtkInformation.h>
#include <vtkTable.h>

#ifdef TTK_ENABLE_TORCH
namespace ttk {
  namespace wnn {

    void makeMatchingVectors(
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &originsMatchings,
      std::vector<mtu::TorchMergeTree<float>> &originsCopy,
      std::vector<mtu::TorchMergeTree<float>> &originsPrimeCopy,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invOriginsMatchingVectorT,
      std::vector<
        std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>>
        &dataMatchings,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
      std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
        &invDataMatchingVectorT,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &reconstMatchings,
      std::vector<std::vector<ttk::ftm::idNode>> &invReconstMatchingVectorT) {
      originsMatchingVectorT.resize(originsMatchings.size());
      invOriginsMatchingVectorT = originsMatchingVectorT;
      for(unsigned int l = 0; l < originsMatchingVectorT.size(); ++l) {
        auto &tree1 = (l == 0 ? originsCopy[0] : originsPrimeCopy[l - 1]);
        auto &tree2 = (l == 0 ? originsPrimeCopy[0] : originsPrimeCopy[l]);
        ttk::axa::getMatchingVector(tree1.mTree, tree2.mTree,
                                    originsMatchings[l],
                                    originsMatchingVectorT[l]);
        ttk::axa::getInverseMatchingVector(tree1.mTree, tree2.mTree,
                                           originsMatchings[l],
                                           invOriginsMatchingVectorT[l]);
      }

      invDataMatchingVectorT.resize(dataMatchings.size());
      for(unsigned int l = 0; l < invDataMatchingVectorT.size(); ++l) {
        invDataMatchingVectorT[l].resize(dataMatchings[l].size());
        for(unsigned int i = 0; i < invDataMatchingVectorT[l].size(); ++i)
          ttk::axa::getInverseMatchingVector(
            originsCopy[l].mTree, recs[i][l].mTree, dataMatchings[l][i],
            invDataMatchingVectorT[l][i]);
      }
      invReconstMatchingVectorT.resize(reconstMatchings.size());
      for(unsigned int i = 0; i < invReconstMatchingVectorT.size(); ++i) {
        auto l = recs[i].size() - 1;
        ttk::axa::getInverseMatchingVector(recs[i][0].mTree, recs[i][l].mTree,
                                           reconstMatchings[i],
                                           invReconstMatchingVectorT[i]);
      }
    }

    void makeDataOutput(
      vtkMultiBlockDataSet *output_data,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
      unsigned int recSize,
      std::vector<vtkDataSet *> &treesSegmentation,
      std::vector<std::vector<double>> &persCorrelationMatrix,
      std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
        &invDataMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invReconstMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<int> &originPersistenceOrder,
      std::vector<vtkUnstructuredGrid *> &treesNodes,
      std::vector<std::vector<int>> &treesNodeCorr,
      std::vector<unsigned int> classId,
      float bestLoss,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      output_data->SetNumberOfBlocks(1);
      vtkSmartPointer<vtkMultiBlockDataSet> data
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      data->SetNumberOfBlocks(1);
      vtkSmartPointer<vtkMultiBlockDataSet> dataSeg
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      dataSeg->SetNumberOfBlocks(recs.size());
      bool outputSegmentation
        = !treesSegmentation.empty() and treesSegmentation[0];
      for(unsigned int l = 0; l < recSize; ++l) {
        vtkSmartPointer<vtkMultiBlockDataSet> out_layer_i
          = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        out_layer_i->SetNumberOfBlocks(recs.size());
        std::vector<ttk::ftm::MergeTree<float> *> trees(recs.size());
        for(unsigned int i = 0; i < recs.size(); ++i)
          trees[i] = &(recs[i][l].mTree);

        // Custom arrays
        std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
          customIntArrays(recs.size());
        std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
          customDoubleArrays(recs.size());
        unsigned int lShift = 0;
        ttk::wnn::computeCustomArrays(
          recs, persCorrelationMatrix, invDataMatchingVectorT,
          invReconstMatchingVectorT, originsMatchingVector,
          originsMatchingVectorT, originsPersPercent, originsPersDiff,
          originPersistenceOrder, l, lShift, customIntArrays,
          customDoubleArrays);
        if(!classId.empty()) {
          for(unsigned int i = 0; i < recs.size(); ++i)
            customIntArrays[i].emplace_back(std::make_tuple(
              "ClassID",
              std::vector<int>(
                recs[i][l].mTree.tree.getNumberOfNodes(), classId[i])));
        }

        // Create output
        if(l == 0)
          ttk::wnn::makeManyOutput(
            trees, treesNodes, treesNodeCorr, out_layer_i, customIntArrays,
            customDoubleArrays, mixtureCoefficient, isPersistenceDiagram,
            convertToDiagram, debugLevel);
        else
          ttk::wnn::makeManyOutput(trees, out_layer_i, customIntArrays,
                                   customDoubleArrays, mixtureCoefficient,
                                   isPersistenceDiagram, convertToDiagram,
                                   debugLevel);
        if(outputSegmentation and l == 0) {
          ttk::wnn::makeManyOutput(
            trees, treesNodes, treesNodeCorr, treesSegmentation, dataSeg,
            customIntArrays, customDoubleArrays, mixtureCoefficient,
            isPersistenceDiagram, convertToDiagram, debugLevel);
        }
        data->SetBlock(l, out_layer_i);
        std::stringstream ss;
        ss << (l == 0 ? "Input" : "Layer") << l;
        data->GetMetaData(l)->Set(vtkCompositeDataSet::NAME(), ss.str());
      }
      output_data->SetBlock(0, data);
      unsigned int num = 0;
      output_data->GetMetaData(num)->Set(
        vtkCompositeDataSet::NAME(), "layersTrees");
      if(outputSegmentation)
        output_data->SetBlock(1, dataSeg);
      vtkNew<vtkFloatArray> lossArray{};
      lossArray->SetName("Loss");
      lossArray->InsertNextTuple1(bestLoss);
      output_data->GetFieldData()->AddArray(lossArray);
    }

    void makeDataOutput(
      vtkMultiBlockDataSet *output_data,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
      unsigned int recSize,
      std::vector<vtkDataSet *> &treesSegmentation,
      std::vector<std::vector<double>> &persCorrelationMatrix,
      std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
        &invDataMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invReconstMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<int> &originPersistenceOrder,
      std::vector<vtkUnstructuredGrid *> &treesNodes,
      std::vector<std::vector<int>> &treesNodeCorr,
      float bestLoss,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      std::vector<unsigned int> classId;
      makeDataOutput(output_data, recs, recSize, treesSegmentation,
                     persCorrelationMatrix, invDataMatchingVectorT,
                     invReconstMatchingVectorT, originsMatchingVectorT,
                     originsMatchingVector, originsPersPercent, originsPersDiff,
                     originPersistenceOrder, treesNodes, treesNodeCorr, classId,
                     bestLoss, mixtureCoefficient, isPersistenceDiagram,
                     convertToDiagram, debugLevel);
    }

    void makeOriginsOutput(
      vtkMultiBlockDataSet *output_origins,
      std::vector<mtu::TorchMergeTree<float>> &originsCopy,
      std::vector<mtu::TorchMergeTree<float>> &originsPrimeCopy,
      std::vector<double> &originPersPercent,
      std::vector<double> &originPersDiff,
      std::vector<int> &originPersistenceOrder,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel) {
      unsigned int noLayers = originsCopy.size();

      output_origins->SetNumberOfBlocks(2);
      // Origins
      vtkSmartPointer<vtkMultiBlockDataSet> origins
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkSmartPointer<vtkMultiBlockDataSet> originsP
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      origins->SetNumberOfBlocks(noLayers);
      originsP->SetNumberOfBlocks(noLayers);
      std::vector<ttk::ftm::MergeTree<float> *> trees(noLayers);
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        customIntArrays(noLayers);
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        customDoubleArrays(noLayers);
      for(unsigned int l = 0; l < noLayers; ++l) {
        trees[l] = &(originsCopy[l].mTree);
        if(l == 0) {
          std::string name2{"OriginPersPercent"};
          customDoubleArrays[l].emplace_back(
            std::make_tuple(name2, originPersPercent));
          std::string name3{"OriginPersDiff"};
          customDoubleArrays[l].emplace_back(
            std::make_tuple(name3, originPersDiff));
          std::string nameOrder{"OriginPersOrder"};
          customIntArrays[l].emplace_back(
            std::make_tuple(nameOrder, originPersistenceOrder));
        }
      }
      ttk::wnn::makeManyOutput(
        trees, origins, customIntArrays, customDoubleArrays, mixtureCoefficient,
        isPersistenceDiagram, convertToDiagram, debugLevel);

      customIntArrays.clear();
      customIntArrays.resize(noLayers);
      customDoubleArrays.clear();
      customDoubleArrays.resize(noLayers);
      for(unsigned int l = 0; l < noLayers; ++l) {
        trees[l] = &(originsPrimeCopy[l].mTree);
        if(l < originsMatchingVector.size()) {
          std::vector<int> customArrayMatching,
            originPersOrder(trees[l]->tree.getNumberOfNodes(), -1);
          for(unsigned int i = 0; i < originsMatchingVector[l].size(); ++i) {
            customArrayMatching.emplace_back(originsMatchingVector[l][i]);
            if(originsMatchingVector[l][i] < originPersistenceOrder.size())
              originPersOrder[i]
                = originPersistenceOrder[originsMatchingVector[l][i]];
          }
          std::string name{"OriginTrueNodeId"};
          customIntArrays[l].emplace_back(
            std::make_tuple(name, customArrayMatching));
          std::string nameOrder{"OriginPersOrder"};
          customIntArrays[l].emplace_back(
            std::make_tuple(nameOrder, originPersOrder));
          std::string name2{"OriginPersPercent"};
          customDoubleArrays[l].emplace_back(
            std::make_tuple(name2, originsPersPercent[l]));
          std::string name3{"OriginPersDiff"};
          customDoubleArrays[l].emplace_back(
            std::make_tuple(name3, originsPersDiff[l]));
        }
      }
      ttk::wnn::makeManyOutput(
        trees, originsP, customIntArrays, customDoubleArrays,
        mixtureCoefficient, isPersistenceDiagram, convertToDiagram, debugLevel);
      output_origins->SetBlock(0, origins);
      output_origins->SetBlock(1, originsP);
      // for(unsigned int l = 0; l < 2; ++l) {
      for(unsigned int l = 0; l < noLayers; ++l) {
        if(l >= 2)
          break;
        std::stringstream ss;
        ss << (l == 0 ? "InputOrigin" : "LayerOrigin") << l;
        auto originsMetaData = origins->GetMetaData(l);
        if(originsMetaData)
          originsMetaData->Set(vtkCompositeDataSet::NAME(), ss.str());
        ss.str("");
        ss << (l == 0 ? "InputOriginPrime" : "LayerOriginPrime") << l;
        auto originsPMetaData = originsP->GetMetaData(l);
        if(originsPMetaData)
          originsPMetaData->Set(vtkCompositeDataSet::NAME(), ss.str());
      }
      unsigned int num = 0;
      output_origins->GetMetaData(num)->Set(
        vtkCompositeDataSet::NAME(), "layersOrigins");
      num = 1;
      output_origins->GetMetaData(num)->Set(
        vtkCompositeDataSet::NAME(), "layersOriginsPrime");
    }

    void makeCoefficientsOutput(
      vtkMultiBlockDataSet *output_coef,
      std::vector<std::vector<torch::Tensor>> &allAlphas,
      std::vector<std::vector<torch::Tensor>> &allScaledAlphas,
      std::vector<std::vector<torch::Tensor>> &allActAlphas,
      std::vector<std::vector<torch::Tensor>> &allActScaledAlphas,
      std::vector<unsigned int> &clusterAsgn,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees) {
      output_coef->SetNumberOfBlocks(allAlphas[0].size());
      for(unsigned int l = 0; l < allAlphas[0].size(); ++l) {
        vtkSmartPointer<vtkTable> coef_table = vtkSmartPointer<vtkTable>::New();
        vtkNew<vtkIntArray> treeIDArray{};
        treeIDArray->SetName("TreeID");
        treeIDArray->SetNumberOfTuples(inputTrees.size());
        for(unsigned int i = 0; i < inputTrees.size(); ++i)
          treeIDArray->SetTuple1(i, i);
        coef_table->AddColumn(treeIDArray);
        auto noVec = allAlphas[0][l].sizes()[0];
        for(unsigned int v = 0; v < noVec; ++v) {
          // Alphas
          vtkNew<vtkFloatArray> tArray{};
          std::string name = ttk::axa::getTableCoefficientName(noVec, v);
          tArray->SetName(name.c_str());
          tArray->SetNumberOfTuples(allAlphas.size());
          // Act Alphas
          vtkNew<vtkFloatArray> actArray{};
          std::string actName = "Act" + name;
          actArray->SetName(actName.c_str());
          actArray->SetNumberOfTuples(allAlphas.size());
          // Scaled Alphas
          vtkNew<vtkFloatArray> tArrayNorm{};
          std::string nameNorm
            = ttk::axa::getTableCoefficientNormName(noVec, v);
          tArrayNorm->SetName(nameNorm.c_str());
          tArrayNorm->SetNumberOfTuples(allAlphas.size());
          // Act Scaled Alphas
          vtkNew<vtkFloatArray> actArrayNorm{};
          std::string actNameNorm = "Act" + nameNorm;
          actArrayNorm->SetName(actNameNorm.c_str());
          actArrayNorm->SetNumberOfTuples(allAlphas.size());
          // Fill Arrays
          for(unsigned int i = 0; i < allAlphas.size(); ++i) {
            tArray->SetTuple1(i, allAlphas[i][l][v].item<float>());
            actArray->SetTuple1(i, allActAlphas[i][l][v].item<float>());
            tArrayNorm->SetTuple1(i, allScaledAlphas[i][l][v].item<float>());
            actArrayNorm->SetTuple1(
              i, allActScaledAlphas[i][l][v].item<float>());
          }
          coef_table->AddColumn(tArray);
          coef_table->AddColumn(actArray);
          coef_table->AddColumn(tArrayNorm);
          coef_table->AddColumn(actArrayNorm);
        }
        if(!clusterAsgn.empty()) {
          vtkNew<vtkIntArray> clusterArray{};
          clusterArray->SetName("ClusterAssignment");
          clusterArray->SetNumberOfTuples(inputTrees.size());
          for(unsigned int i = 0; i < clusterAsgn.size(); ++i)
            clusterArray->SetTuple1(i, clusterAsgn[i]);
          coef_table->AddColumn(clusterArray);
        }
        if(l == 0) {
          vtkNew<vtkIntArray> treesNoNodesArray{};
          treesNoNodesArray->SetNumberOfTuples(recs.size());
          treesNoNodesArray->SetName("treeNoNodes");
          for(unsigned int i = 0; i < recs.size(); ++i)
            treesNoNodesArray->SetTuple1(
              i, recs[i][0].mTree.tree.getNumberOfNodes());
          coef_table->AddColumn(treesNoNodesArray);
        }
        output_coef->SetBlock(l, coef_table);
        std::stringstream ss;
        ss << "Coef" << l;
        output_coef->GetMetaData(l)->Set(vtkCompositeDataSet::NAME(), ss.str());
      }

      // Copy Field Data
      // - aggregate input field data
      for(unsigned int b = 0; b < inputTrees[0]->GetNumberOfBlocks(); ++b) {
        vtkNew<vtkFieldData> fd{};
        fd->CopyStructure(inputTrees[0]->GetBlock(b)->GetFieldData());
        fd->SetNumberOfTuples(inputTrees.size());
        for(size_t i = 0; i < inputTrees.size(); ++i) {
          fd->SetTuple(i, 0, inputTrees[i]->GetBlock(b)->GetFieldData());
        }

        // - copy input field data to output row data
        for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
          auto array = fd->GetAbstractArray(i);
          array->SetName(array->GetName());
          vtkTable::SafeDownCast(output_coef->GetBlock(0))->AddColumn(array);
        }
      }
    }

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

  } // namespace wnn
} // namespace ttk
#endif
