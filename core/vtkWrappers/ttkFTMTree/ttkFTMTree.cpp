#include <vtkCellData.h>
#include <ttkFTMTree.h>

vtkStandardNewMacro(ttkFTMTree)

    ttkFTMTree::ttkFTMTree()
    : ScalarField{},
      UseInputOffsetScalarField{},
      InputOffsetScalarFieldName{},
      ScalarFieldId{},
      OffsetFieldId{-1},

      triangulation_{},
      inputScalars_{},
      offsets_{},
      inputOffsets_{},
      hasUpdatedMesh_{}
{
   SetNumberOfInputPorts(1);
   SetNumberOfOutputPorts(3);
}

ttkFTMTree::~ttkFTMTree()
{
   if (offsets_)
      offsets_->Delete();
}

int ttkFTMTree::FillInputPortInformation(int port, vtkInformation* info)
{
   if (port == 0)
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
   return 1;
}

int ttkFTMTree::FillOutputPortInformation(int port, vtkInformation* info)
{
   switch (port) {
      case 0:
      case 1:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
         break;

      case 2:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
         break;
   }

   return 1;
}

int ttkFTMTree::setupTriangulation(vtkDataSet* input)
{
   triangulation_ = ttkTriangulation::getTriangulation(input);
#ifndef withKamikaze
   if (!triangulation_) {
      cerr << "[ttkFTMTree] Error : ttkTriangulation::getTriangulation() is null." << endl;
      return -1;
   }
#endif

   hasUpdatedMesh_ = ttkTriangulation::hasChangedConnectivity(triangulation_, input, this);

   triangulation_->setWrapper(this);
   contourForests_.setDebugLevel(debugLevel_);
   contourForests_.setThreadNumber(threadNumber_);
   contourForests_.setupTriangulation(triangulation_);

#ifndef withKamikaze
   if (triangulation_->isEmpty()) {
      cerr << "[ttkFTMTree] Error : ttkTriangulation allocation problem." << endl;
      return -1;
   }
#endif

   return 0;
}

int ttkFTMTree::getScalars(vtkDataSet* input)
{
   vtkPointData* pointData = input->GetPointData();

#ifndef withKamikaze
   if (!pointData) {
      cerr << "[ttkFTMTree] Error : input has no point data." << endl;
      return -1;
   }
#endif

   if (ScalarField.length()) {
      inputScalars_ = pointData->GetArray(ScalarField.data());
   } else {
      inputScalars_ = pointData->GetArray(ScalarFieldId);
      if (inputScalars_)
         ScalarField = inputScalars_->GetName();
   }

#ifndef withKamikaze
   if (!inputScalars_) {
      cerr << "[ttkFTMTree] Error : input scalar field pointer is null." << endl;
      return -3;
   }
#endif

   return 0;
}

int ttkFTMTree::getOffsets(vtkDataSet* input)
{
   if (OffsetFieldId != -1) {
      inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
      if (inputOffsets_) {
         InputOffsetScalarFieldName = inputOffsets_->GetName();
         UseInputOffsetScalarField  = true;
      }
   }

   if (UseInputOffsetScalarField and InputOffsetScalarFieldName.length())
      inputOffsets_ = input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
   else {
      if (hasUpdatedMesh_ and offsets_) {
         offsets_->Delete();
         offsets_ = nullptr;
      }

      if (!offsets_) {
         const int numberOfVertices = input->GetNumberOfPoints();

         offsets_ = vtkIntArray::New();
         offsets_->SetNumberOfComponents(1);
         offsets_->SetNumberOfTuples(numberOfVertices);
         offsets_->SetName("OffsetsScalarField");
         for (int i = 0; i < numberOfVertices; ++i)
            offsets_->SetTuple1(i, i);
      }

      inputOffsets_ = offsets_;
   }

#ifndef withKamikaze
   if (!inputOffsets_) {
      cerr << "[ttkFTMTree] Error : wrong input offset scalar field." << endl;
      return -1;
   }
#endif

   return 0;
}

NodeType ttkFTMTree::getNodeType(const Node* node)
{
   int upDegree{};
   int downDegree{};
   if (treeType_ == TreeType::Join or treeType_ == TreeType::Contour) {
      upDegree   = node->getNumberOfUpSuperArcs();
      downDegree = node->getNumberOfDownSuperArcs();
   } else {
      downDegree = node->getNumberOfUpSuperArcs();
      upDegree   = node->getNumberOfDownSuperArcs();
   }
   int degree = upDegree + downDegree;

   // saddle point
   if (degree > 1) {
      if (upDegree == 2 and downDegree == 1)
         return NodeType::Saddle2;
      else if (upDegree == 1 && downDegree == 2)
         return NodeType::Saddle1;
      else if (upDegree == 1 && downDegree == 1)
         return NodeType::Regular;
      else
         return NodeType::Degenerate;
   }
   // local extremum
   else {
      if (upDegree)
         return NodeType::Local_minimum;
      else
         return NodeType::Local_maximum;
   }
}

int ttkFTMTree::getSkeletonNodes(FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonNodes)
{
   vtkSmartPointer<vtkUnstructuredGrid> skeletonNodes = vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifndef withKamikaze
   if (!skeletonNodes) {
      cerr << "[ttkFTMTree] Error : vtkUnstructuredGrid allocation problem." << endl;
      return -1;
   }
#endif

   const idVertex numberOfNodes = tree->getNumberOfNodes();
#ifndef withKamikaze
   if (!numberOfNodes) {
      cerr << "[ttkFTMTree] Error : tree has no nodes." << endl;
      return -2;
   }
#endif

   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
   if (!points) {
      cerr << "[ttkFTMTree] Error : vtkPoints allocation problem." << endl;
      return -3;
   }
#endif

   vtkSmartPointer<vtkIntArray> nodeIds = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!nodeIds) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem." << endl;
      return -4;
   }
#endif
   nodeIds->SetNumberOfComponents(1);
   nodeIds->SetName("NodeId");

   vtkSmartPointer<vtkIntArray> vertexIds = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!vertexIds) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem" << endl;
      return -5;
   }
#endif
   vertexIds->SetNumberOfComponents(1);
   vertexIds->SetName("VertexId");

   vtkSmartPointer<vtkIntArray> nodeTypes = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!nodeTypes) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem" << endl;
      return -6;
   }
#endif
   nodeTypes->SetNumberOfComponents(1);
   nodeTypes->SetName("NodeType");

   for (idVertex nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
      const Node* node = tree->getNode(nodeId);
#ifndef withKamikaze
      if (!node) {
         cerr << "[ttkFTMTree] Error : node " << nodeId << " is null." << endl;
         return -7;
      }
#endif

      const idVertex vertexId = node->getVertexId();
      const idVertex nodeType = static_cast<idVertex>(getNodeType(node));

      float point[3];
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);

      points->InsertNextPoint(point);
      nodeIds->InsertNextTuple1(nodeId);
      vertexIds->InsertNextTuple1(vertexId);
      nodeTypes->InsertNextTuple1(nodeType);
   }

   skeletonNodes->SetPoints(points);

   vtkPointData* pointData = skeletonNodes->GetPointData();
#ifndef withKamikaze
   if (!pointData) {
      cerr << "[ttkFTMTree] Error : output skeleton nodes has no point data." << endl;
      return -8;
   }
#endif

   pointData->AddArray(nodeIds);
   pointData->AddArray(vertexIds);
   pointData->AddArray(nodeTypes);

   outputSkeletonNodes->ShallowCopy(skeletonNodes);

   return 0;
}

int ttkFTMTree::addDirectSkeletonArc(FTMTree_MT* tree, SuperArc* arc, vtkPoints* points,
                                        vtkUnstructuredGrid* skeletonArcs)
{
   float     point[3];
   vtkIdType ids[2];

   const idVertex upNodeId   = arc->getUpNodeId();
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);
   ids[0] = points->InsertNextPoint(point);

   const idVertex downNodeId   = arc->getDownNodeId();
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);
   ids[1] = points->InsertNextPoint(point);

   skeletonArcs->InsertNextCell(VTK_LINE, 2, ids);

   return 0;
}

int ttkFTMTree::addSampledSkeletonArc(FTMTree_MT* tree, SuperArc* arc, const int samplingLevel,
                                         vtkPoints* points, vtkUnstructuredGrid* skeletonArcs)
{
   float     point[3];
   vtkIdType ids[2];

   const idVertex downNodeId   = tree->getLowerNodeId(arc);
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);
   const vtkIdType downId    = points->InsertNextPoint(point);
   const double    scalarMin = inputScalars_->GetTuple1(downVertexId);

   const idVertex upNodeId   = tree->getUpperNodeId(arc);
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);
   const vtkIdType upId      = points->InsertNextPoint(point);
   const double    scalarMax = inputScalars_->GetTuple1(upVertexId);

   const double delta       = (scalarMax - scalarMin) / (samplingLevel + 1);
   double       scalarLimit = scalarMin + delta;

   ids[0] = downId;

   int       c = 0;
   float     sum[3]{0, 0, 0};
   for (const idVertex vertexId : *arc) {
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);
      const double scalarVertex = inputScalars_->GetTuple1(vertexId);

      if (scalarVertex < scalarLimit) {
         sum[0] += point[0];
         sum[1] += point[1];
         sum[2] += point[2];
         ++c;
      } else {
         sum[0] /= c;
         sum[1] /= c;
         sum[2] /= c;

         ids[1] = points->InsertNextPoint(sum);
         skeletonArcs->InsertNextCell(VTK_LINE, 2, ids);
         ids[0] = ids[1];

         scalarLimit += delta;
         sum[0] = 0;
         sum[1] = 0;
         sum[2] = 0;
         c      = 0;
      }
   }

   ids[1] = upId;

   skeletonArcs->InsertNextCell(VTK_LINE, 2, ids);

   return 0;
}

int ttkFTMTree::addCompleteSkeletonArc(FTMTree_MT* tree, SuperArc* arc, vtkPoints* points,
                                          vtkUnstructuredGrid* skeletonArcs)
{
   float     point[3];
   vtkIdType ids[2];

   const idVertex downNodeId   = arc->getDownNodeId();
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);
   ids[0] = points->InsertNextPoint(point);

   for (const idVertex vertexId : *arc) {
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);
      ids[1] = points->InsertNextPoint(point);

      skeletonArcs->InsertNextCell(VTK_LINE, 2, ids);

      ids[0] = ids[1];
   }

   const idVertex upNodeId   = arc->getUpNodeId();
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);
   ids[1] = points->InsertNextPoint(point);

   skeletonArcs->InsertNextCell(VTK_LINE, 2, ids);

   return 0;
}

int ttkFTMTree::getSkeletonArcs(FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonArcs)
{
   vtkSmartPointer<vtkUnstructuredGrid> skeletonArcs = vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifndef withKamikaze
   if (!skeletonArcs) {
      cerr << "[ttkFTMTree] Error : vtkUnstructuredGrid allocation problem." << endl;
      return -1;
   }
#endif

   const idVertex numberOfSuperArcs = tree->getNumberOfSuperArcs();
#ifndef withKamikaze
   if (!numberOfSuperArcs) {
      cerr << "[ttkFTMTree] Error : tree has no super arcs." << endl;
      return -2;
   }
#endif

#ifdef withStatsTime
   vtkSmartPointer<vtkFloatArray> startArcs=vtkSmartPointer<vtkFloatArray>::New();
   startArcs->SetName("Start");
# ifndef withKamikaze
   if (!startArcs) {
      cerr << "[ttkFTMTree] Error : vtkFloatArray start allocation problem." << endl;
      return -5;
   }
# endif

   vtkSmartPointer<vtkFloatArray> endArcs=vtkSmartPointer<vtkFloatArray>::New();
   endArcs->SetName("End");
# ifndef withKamikaze
   if (!endArcs) {
      cerr << "[ttkFTMTree] Error : vtkFloatArray end allocation problem." << endl;
      return -5;
   }
# endif

   vtkSmartPointer<vtkFloatArray> timeArcs=vtkSmartPointer<vtkFloatArray>::New();
   timeArcs->SetName("Time");
# ifndef withKamikaze
   if (!timeArcs) {
      cerr << "[ttkFTMTree] Error : vtkFloatArray time allocation problem." << endl;
      return -5;
   }
# endif

   vtkSmartPointer<vtkIntArray> origArcs=vtkSmartPointer<vtkIntArray>::New();
   origArcs->SetName("Origin");
# ifndef withKamikaze
   if (!origArcs) {
      cerr << "[ttkFTMTree] Error : vtkIntArray origin allocation problem." << endl;
      return -6;
   }
# endif

   vtkSmartPointer<vtkIntArray> tasksArcs=vtkSmartPointer<vtkIntArray>::New();
   tasksArcs->SetName("Tasks");
# ifndef withKamikaze
   if (!tasksArcs) {
      cerr << "[ttkFTMTree] Error : vtkIntArray tasks allocation problem." << endl;
      return -5;
   }
# endif

#endif

   vtkSmartPointer<vtkIntArray> sizeArcs=vtkSmartPointer<vtkIntArray>::New();
   sizeArcs->SetName("Size");
#ifndef withKamikaze
   if(!sizeArcs) {
      cerr << "[ttkFTMTree] Error : vtkIntArray size allocation problem." << endl;
      return -7;
   }
#endif

   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
   if (!points) {
      cerr << "[ttkFTMTree] Error : vtkPoints allocation problem." << endl;
      return -8;
   }
#endif

   const int samplingLevel = SuperArcSamplingLevel;

   for (idVertex i = 0; i < numberOfSuperArcs; ++i) {
      SuperArc* arc = tree->getSuperArc(i);

      const int numberOfRegularNodes = arc->getNumberOfRegularNodes();
      if (numberOfRegularNodes > 0 and samplingLevel > 0)
          addSampledSkeletonArc(tree, arc, samplingLevel, points, skeletonArcs);
      else if (samplingLevel == -1)
          addCompleteSkeletonArc(tree, arc, points, skeletonArcs);
      else
          addDirectSkeletonArc(tree, arc, points, skeletonArcs);

#ifdef withStatsTime
      startArcs->InsertNextTuple1(tree->getArcStart(i));
      endArcs->InsertNextTuple1(tree->getArcEnd(i));
      timeArcs->InsertNextTuple1(tree->getArcEnd(i) - tree->getArcStart(i));
      origArcs->InsertNextTuple1(tree->getArcOrig(i));
      tasksArcs->InsertNextTuple1(tree->getArcActiveTasks(i));
#endif

      sizeArcs->InsertNextTuple1(tree->getArcSize(i));
   }

   skeletonArcs->SetPoints(points);
#ifdef withStatsTime
   skeletonArcs->GetCellData()->AddArray(startArcs);
   skeletonArcs->GetCellData()->AddArray(endArcs);
   skeletonArcs->GetCellData()->AddArray(timeArcs);
   skeletonArcs->GetCellData()->AddArray(origArcs);
   skeletonArcs->GetCellData()->AddArray(tasksArcs);
#endif
   skeletonArcs->GetCellData()->AddArray(sizeArcs);
   outputSkeletonArcs->ShallowCopy(skeletonArcs);

   return 0;
}

int ttkFTMTree::getSegmentation(FTMTree_MT* tree, vtkDataSet* input,
                                   vtkDataSet* outputSegmentation)
{
   outputSegmentation->ShallowCopy(input);

   const idVertex numberOfSuperArcs = tree->getNumberOfSuperArcs();
#ifndef withKamikaze
   if (!numberOfSuperArcs) {
      cerr << "[ttkFTMTree] Error : tree has no super arcs." << endl;
      return -1;
   }
#endif

   const idVertex numberOfVertices = triangulation_->getNumberOfVertices();
#ifndef withKamikaze
   if (!numberOfVertices) {
      cerr << "[ttkFTMTree] Error : triangulation has no vertices." << endl;
      return -2;
   }
#endif

   // field
   vtkSmartPointer<vtkIntArray> regionIds = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!regionIds) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem." << endl;
      return -3;
   }
#endif
   regionIds->SetName("SegmentationId");
   regionIds->SetNumberOfComponents(1);
   regionIds->SetNumberOfTuples(numberOfVertices);

   vtkSmartPointer<vtkIntArray> regionTypes = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!regionTypes) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem." << endl;
      return -4;
   }
#endif
   regionTypes->SetName("RegionType");
   regionTypes->SetNumberOfComponents(1);
   regionTypes->SetNumberOfTuples(numberOfVertices);

   vtkSmartPointer<vtkIntArray> regionSizes = vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
   if (!regionSizes) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem." << endl;
      return -5;
   }
#endif
   regionSizes->SetName("RegionSize");
   regionSizes->SetNumberOfComponents(1);
   regionSizes->SetNumberOfTuples(numberOfVertices);

   vtkSmartPointer<vtkDoubleArray> regionSpans = vtkSmartPointer<vtkDoubleArray>::New();
#ifndef withKamikaze
   if (!regionSpans) {
      cerr << "[ttkFTMTree] Error : vtkIntArray allocation problem." << endl;
      return -6;
   }
#endif
   regionSpans->SetName("RegionSpan");
   regionSpans->SetNumberOfComponents(1);
   regionSpans->SetNumberOfTuples(numberOfVertices);

   int id{};

   // arcs
   for (idVertex arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
      SuperArc* arc = tree->getSuperArc(arcId);
#ifndef withKamikaze
      if (!arc) {
         cerr << "[ttkFTMTree] Error : FTMTree arc " << arcId << " is null." << endl;
         return -7;
      }
#endif

      const int   upNodeId = arc->getUpNodeId();
      const Node* upNode   = tree->getNode(upNodeId);
#ifndef withKamikaze
      if (!upNode) {
          cerr << "[ttkFTMTree] Error : FTMTree node" << upNodeId << " is null." << endl;
          return -8;
      }
#endif
      const NodeType upNodeType = getNodeType(upNode);
      const int      upVertex   = upNode->getVertexId();
      float          coordUp[3];
      triangulation_->getVertexPoint(upVertex, coordUp[0], coordUp[1], coordUp[2]);

      const int   downNodeId = arc->getDownNodeId();
      const Node* downNode   = tree->getNode(downNodeId);
#ifndef withKamikaze
      if (!downNode) {
          cerr << "[ttkFTMTree] Error : FTMTree node" << downNodeId << " is null." << endl;
          return -9;
      }
#endif
      const NodeType downNodeType = getNodeType(downNode);
      const int      downVertex   = downNode->getVertexId();
      float          coordDown[3];
      triangulation_->getVertexPoint(downVertex, coordDown[0], coordDown[1], coordDown[2]);

      const int    regionSize = tree->getSuperArc(arcId)->getNumberOfRegularNodes();
      const double regionSpan = Geometry::distance(coordUp, coordDown);

      int regionType{};
      if (upNodeType == NodeType::Local_minimum or downNodeType == NodeType::Local_minimum)
          regionType = static_cast<int>(ArcType::Min_arc);
      else if (upNodeType == NodeType::Local_maximum or downNodeType == NodeType::Local_maximum)
          regionType = static_cast<int>(ArcType::Max_arc);
      else if (upNodeType == NodeType::Saddle1 and downNodeType == NodeType::Saddle1)
          regionType = static_cast<int>(ArcType::Saddle1_arc);
      else if (upNodeType == NodeType::Saddle2 and downNodeType == NodeType::Saddle2)
          regionType = static_cast<int>(ArcType::Saddle2_arc);
      else
          regionType = static_cast<int>(ArcType::Saddle1_saddle2_arc);

      // critical points
      regionIds->SetTuple1(upVertex, id);
      regionIds->SetTuple1(downVertex, id);
      regionSizes->SetTuple1(upVertex, regionSize);
      regionSizes->SetTuple1(downVertex, regionSize);
      regionSpans->SetTuple1(upVertex, regionSpan);
      regionSpans->SetTuple1(downVertex, regionSpan);
      regionTypes->SetTuple1(upVertex, -1);
      regionTypes->SetTuple1(downVertex, -1);

      // regular nodes
      for (const idVertex vertexId : *arc) {
          regionIds->SetTuple1(vertexId, id);
          regionSizes->SetTuple1(vertexId, regionSize);
          regionSpans->SetTuple1(vertexId, regionSpan);
          regionTypes->SetTuple1(vertexId, regionType);
      }

      ++id;
   }

   vtkPointData* pointData = outputSegmentation->GetPointData();
#ifndef withKamikaze
   if (!pointData) {
      cerr << "[ttkFTMTree] Error : output segmentation has no point data." << endl;
      return -9;
   }
#endif

   pointData->AddArray(regionIds);
   pointData->AddArray(regionSizes);
   pointData->AddArray(regionSpans);
   pointData->AddArray(regionTypes);

   return 0;
}

int ttkFTMTree::doIt(vector<vtkDataSet*>& inputs, vector<vtkDataSet*>& outputs)
{
   Memory m;

   vtkDataSet*          input               = inputs[0];
   vtkUnstructuredGrid* outputSkeletonNodes = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
   vtkUnstructuredGrid* outputSkeletonArcs  = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
   vtkDataSet*          outputSegmentation  = outputs[2];

   if (setupTriangulation(input)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong triangulation." << endl;
      return -1;
#endif
   }

   const idVertex numberOfVertices = triangulation_->getNumberOfVertices();
#ifndef withKamikaze
   if (!numberOfVertices) {
      cerr << "[ttkFTMTree] Error : input data has no vertices." << endl;
      return -2;
   }
#endif

   if (getScalars(input)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong scalars." << endl;
      return -3;
#endif
   }

   if (getOffsets(input)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong offsets." << endl;
      return -4;
#endif
   }

   if(debugLevel_) {
       cout << "Launch on field : " << ScalarField << endl;
   }

   vector<idVertex> offsets(numberOfVertices);
   for (idVertex i = 0; i < numberOfVertices; ++i)
      offsets[i]   = inputOffsets_->GetTuple1(i);

   contourForests_.setVertexScalars(inputScalars_->GetVoidPointer(0));
   // for now, gets only vector<idVertex>, no raw data
   // contour.setVertexSoSoffsets(inputOffsets_->GetVoidPointer(0));
   contourForests_.setVertexSoSoffsets(offsets);
   contourForests_.setTreeType(treeType_);

   switch (inputScalars_->GetDataType()) {
      vtkTemplateMacro(({ contourForests_.build<VTK_TT>(); }));
   }

   FTMTree_MT* tree{};
   switch (treeType_) {
      case TreeType::Join:
         tree = contourForests_.getJoinTree();
         break;
      case TreeType::Split:
         tree = contourForests_.getSplitTree();
         break;
      case TreeType::Contour:
         tree = &contourForests_;
         break;
   }
#ifndef withKamikaze
   if (!tree) {
      cerr << "[ttkFTMTree] tree is null." << endl;
      return -6;
   }
#endif

   if (getSkeletonNodes(tree, outputSkeletonNodes)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong properties on skeleton nodes." << endl;
      return -7;
#endif
   }

   if (getSkeletonArcs(tree, outputSkeletonArcs)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong properties on skeleton arcs." << endl;
      return -8;
#endif
   }

   if (getSegmentation(tree, input, outputSegmentation)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong properties on segmentation." << endl;
      return -9;
#endif
   }

   {
      stringstream msg;
      msg << "[ttkFTMTree] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}
