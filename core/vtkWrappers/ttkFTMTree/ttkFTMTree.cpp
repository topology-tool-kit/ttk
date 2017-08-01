#include <ttkFTMTree.h>

vtkStandardNewMacro(ttkFTMTree);

ttkFTMTree::ttkFTMTree()
    : ScalarField{},
      UseInputOffsetScalarField{},
      InputOffsetScalarFieldName{},
      ScalarFieldId{},
      OffsetFieldId{-1},
      params_{},
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

   triangulation_->setWrapper(this);
   ftmTree_.setDebugLevel(debugLevel_);
   ftmTree_.setThreadNumber(threadNumber_);
   ftmTree_.setupTriangulation(triangulation_);

   hasUpdatedMesh_ = ttkTriangulation::hasChangedConnectivity(triangulation_, input, this);

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
   if (GetTreeType() == TreeType::Join or GetTreeType() == TreeType::Contour) {
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

int ttkFTMTree::addDirectSkeletonArc(FTMTree_MT* tree, idSuperArc arcId, vtkPoints* points,
                                     vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   SuperArc* arc = tree->getSuperArc(arcId);
   float     point[3];
   vtkIdType pointIds[2];

   const idVertex downNodeId   = tree->getLowerNodeId(arc);
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);
   // Get or create first point of the arc
   if (arcData.pointIds[downNodeId] == nullNodes) {
      const vtkIdType downId       = points->InsertNextPoint(point);
      pointIds[0]                  = downId;
      arcData.pointIds[downNodeId] = downId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      pointIds[0] = arcData.pointIds[downNodeId];
   }

   const idVertex upNodeId   = tree->getUpperNodeId(arc);
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);
   // Get or create last point of the arc
   if (arcData.pointIds[upNodeId] == nullNodes) {
      const vtkIdType upId       = points->InsertNextPoint(point);
      pointIds[1]                = upId;
      arcData.pointIds[upNodeId] = upId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      pointIds[1] = arcData.pointIds[upNodeId];
   }

   skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

   arcData.fillArray(arcId, triangulation_, tree, params_);

   return 0;
}

int ttkFTMTree::addSampledSkeletonArc(FTMTree_MT* tree, idSuperArc arcId, const int samplingLevel,
                                      vtkPoints* points, vtkUnstructuredGrid* skeletonArcs,
                                      ArcData& arcData)
{
   SuperArc* arc = tree->getSuperArc(arcId);
   float     point[3];
   vtkIdType pointIds[2];

   const idVertex downNodeId   = tree->getLowerNodeId(arc);
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);
   const double scalarMin = inputScalars_->GetTuple1(downVertexId);

   // Get or create first point of the arc
   vtkIdType downId;
   if (arcData.pointIds[downNodeId] == nullNodes) {
      downId                       = points->InsertNextPoint(point);
      arcData.pointIds[downNodeId] = downId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      downId = arcData.pointIds[downNodeId];
   }

   const idVertex upNodeId   = tree->getUpperNodeId(arc);
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);
   const double    scalarMax = inputScalars_->GetTuple1(upVertexId);

   const double delta       = (scalarMax - scalarMin) / (samplingLevel + 1);
   double       scalarLimit = scalarMin + delta;

   // Get or create last point of the arc
   vtkIdType upId;
   if (arcData.pointIds[upNodeId] == nullNodes) {
      upId                       = points->InsertNextPoint(point);
      arcData.pointIds[upNodeId] = upId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      upId = arcData.pointIds[upNodeId];
   }

   pointIds[0] = downId;

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
         if (c) {
            sum[0] /= c;
            sum[1] /= c;
            sum[2] /= c;

            pointIds[1] = points->InsertNextPoint(sum);
            arcData.regularMask->InsertNextTuple1(true);
            skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

            arcData.fillArray(arcId, triangulation_, tree, params_);

            pointIds[0] = pointIds[1];
         }

         scalarLimit += delta;
         sum[0] = 0;
         sum[1] = 0;
         sum[2] = 0;
         c      = 0;
      }
   }

   pointIds[1] = upId;

   skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

   arcData.fillArray(arcId, triangulation_, tree, params_);

   return 0;
}

int ttkFTMTree::addCompleteSkeletonArc(FTMTree_MT* tree, idSuperArc arcId, vtkPoints* points,
                                       vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   SuperArc* arc = tree->getSuperArc(arcId);
   float     point[3];
   vtkIdType pointIds[2];

   const idVertex downNodeId   = tree->getLowerNodeId(arc);
   const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
   triangulation_->getVertexPoint(downVertexId, point[0], point[1], point[2]);

   // Get or create first point of the arc
   vtkIdType downId;
   if (arcData.pointIds[downNodeId] == nullNodes) {
      downId                       = points->InsertNextPoint(point);
      arcData.pointIds[downNodeId] = downId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      downId = arcData.pointIds[downNodeId];
   }

   pointIds[0] = downId;

   for (const idVertex vertexId : *arc) {
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);
      pointIds[1] = points->InsertNextPoint(point);
      arcData.regularMask->InsertNextTuple1(true);

      skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

      arcData.fillArray(arcId, triangulation_, tree, params_);

      pointIds[0] = pointIds[1];
   }

   const idVertex upNodeId   = tree->getUpperNodeId(arc);
   const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
   triangulation_->getVertexPoint(upVertexId, point[0], point[1], point[2]);

   // Get or create last point of the arc
   vtkIdType upId;
   if (arcData.pointIds[upNodeId] == nullNodes) {
      upId                       = points->InsertNextPoint(point);
      arcData.pointIds[upNodeId] = upId;
      arcData.regularMask->InsertNextTuple1(false);
   } else {
      upId = arcData.pointIds[upNodeId];
   }

   pointIds[1] = upId;

   skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

   arcData.fillArray(arcId, triangulation_, tree, params_);

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

   ArcData arcData;
   arcData.init(tree->getNumberOfNodes(), params_);

   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
   if (!points) {
      cerr << "[ttkFTMTree] Error : vtkPoints allocation problem." << endl;
      return -8;
   }
#endif

   const int samplingLevel = SuperArcSamplingLevel;

   for (idVertex arcId = 0; arcId < numberOfSuperArcs; ++arcId) {

      const int numberOfRegularNodes = tree->getArcSize(arcId);
      if (numberOfRegularNodes > 0 and samplingLevel > 0) {
         addSampledSkeletonArc(tree, arcId, samplingLevel, points, skeletonArcs, arcData);
      } else if (samplingLevel == -1) {
         addCompleteSkeletonArc(tree, arcId, points, skeletonArcs, arcData);
      } else {
         addDirectSkeletonArc(tree, arcId, points, skeletonArcs, arcData);
      }
   }

   skeletonArcs->SetPoints(points);
   arcData.addArray(skeletonArcs, params_);
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

   VertData vertdata;

   vertdata.init(numberOfVertices, params_);

   // arcs
   for (idVertex arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
       vertdata.fillArray(arcId, triangulation_, tree, params_);
   }

   vtkPointData* pointData = outputSegmentation->GetPointData();
#ifndef withKamikaze
   if (!pointData) {
      cerr << "[ttkFTMTree] Error : output segmentation has no point data." << endl;
      return -9;
   }
#endif

   vertdata.addArray(pointData, params_);

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

   ftmTree_.setVertexScalars(inputScalars_->GetVoidPointer(0));
   ftmTree_.setVertexSoSoffsets(offsets);
   ftmTree_.setTreeType(GetTreeType());
   ftmTree_.setSegmentation(GetWithSegmentation());
   ftmTree_.setNormalizeIds(GetWithNormalize());

   switch (inputScalars_->GetDataType()) {
      vtkTemplateMacro(({ ftmTree_.build<VTK_TT>(); }));
   }

   UpdateProgress(0.50);

   FTMTree_MT* tree = ftmTree_.getTree(GetTreeType());
#ifndef withKamikaze
   if (!tree) {
      cerr << "[ttkFTMTree] tree is null." << endl;
      return -6;
   }
#endif

   UpdateProgress(0.70);

   if (getSkeletonNodes(tree, outputSkeletonNodes)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong properties on skeleton nodes." << endl;
      return -7;
#endif
   }

   UpdateProgress(0.75);

   if (getSkeletonArcs(tree, outputSkeletonArcs)) {
#ifndef withKamikaze
      cerr << "[ttkFTMTree] Error : wrong properties on skeleton arcs." << endl;
      return -8;
#endif
   }

   if (GetWithSegmentation()) {
      if (getSegmentation(tree, input, outputSegmentation)) {
#ifndef withKamikaze
         cerr << "[ttkFTMTree] Error : wrong properties on segmentation." << endl;
         return -9;
#endif
      }
   }

   UpdateProgress(1);

   {
      stringstream msg;
      msg << "[ttkFTMTree] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}
