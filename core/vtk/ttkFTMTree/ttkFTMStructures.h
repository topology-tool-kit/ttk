#ifndef TTKFTMSTRUCTURES_H
#define TTKFTMSTRUCTURES_H

#include <FTMTree.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkCharArray.h>

struct LocalFTM {
   ttk::ftm::FTMTree tree;
   ttk::ftm::idNode  offset;
};

struct WrapperData {
   template <typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> initArray(const char* fieldName, size_t nbElmnt)
   {
      vtkSmartPointer<vtkArrayType> arr = vtkSmartPointer<vtkArrayType>::New();
      arr->SetName(fieldName);
      arr->SetNumberOfComponents(1);
      arr->SetNumberOfTuples(nbElmnt);

#ifndef TTK_ENABLE_KAMIKAZE
      if (!arr) {
         std::cerr << "[ttkFTMTree] Error, unable to allocate " << fieldName
              << " the program will likely crash" << std::endl;
      }
#endif
      return arr;
   }

   inline static ttk::ftm::NodeType getNodeType(ttk::ftm::FTMTree_MT&  tree,
                                                const ttk::ftm::idNode nodeId,
                                                ttk::ftm::Params       params)
   {
      const ttk::ftm::Node* node = tree.getNode(nodeId);
      int upDegree{};
      int downDegree{};
      if (params.treeType == ttk::ftm::TreeType::Join or
          params.treeType == ttk::ftm::TreeType::Contour) {
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
            return ttk::ftm::NodeType::Saddle2;
         else if (upDegree == 1 and downDegree == 2)
            return ttk::ftm::NodeType::Saddle1;
         else if (upDegree == 1 and downDegree == 1)
            return ttk::ftm::NodeType::Regular;
         else
            return ttk::ftm::NodeType::Degenerate;
      }
      // local extremum
      else {
         if (upDegree)
            return ttk::ftm::NodeType::Local_minimum;
         else
            return ttk::ftm::NodeType::Local_maximum;
      }
   }
};

struct ArcData : public WrapperData {
   std::vector<ttk::SimplexId>          point_ids;
   vtkSmartPointer<vtkCharArray>   point_regularMask;
   vtkSmartPointer<vtkFloatArray>  point_scalar;
#ifdef TTK_USE_64BIT_IDS
   vtkSmartPointer<vtkIdTypeArray>    cell_ids;
   vtkSmartPointer<vtkIdTypeArray>    cell_upNode;
   vtkSmartPointer<vtkIdTypeArray>    cell_downNode;
   vtkSmartPointer<vtkIdTypeArray>    cell_sizeArcs;
#else
   vtkSmartPointer<vtkIntArray>    cell_ids;
   vtkSmartPointer<vtkIntArray>    cell_upNode;
   vtkSmartPointer<vtkIntArray>    cell_downNode;
   vtkSmartPointer<vtkIntArray>    cell_sizeArcs;
#endif
   vtkSmartPointer<vtkDoubleArray> cell_spanArcs;

   inline int init(std::vector<LocalFTM>& ftmTree, ttk::ftm::Params params)
   {
      ttk::ftm::idSuperArc nbArcs       = 0;
      ttk::ftm::idSuperArc nbNodes      = 0;
      ttk::ftm::idSuperArc samplePoints = 0;
      ttk::ftm::idVertex   nbVerts      = 0;

      for (auto& t : ftmTree) {
         ttk::ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         nbArcs += tree->getNumberOfSuperArcs();
         nbNodes += tree->getNumberOfNodes();
         samplePoints += params.samplingLvl >= 0
                             ? tree->getNumberOfNodes() + (nbArcs * params.samplingLvl)
                             : tree->getNumberOfVertices();
         nbVerts += tree->getNumberOfVertices();
      }

      point_ids.resize(nbVerts, ttk::ftm::nullVertex);
#ifdef TTK_USE_64BIT_IDS
      cell_ids          = initArray<vtkIdTypeArray>("SegmentationId", samplePoints);
      cell_upNode       = initArray<vtkIdTypeArray>("upNodeId", samplePoints);
      cell_downNode     = initArray<vtkIdTypeArray>("downNodeId", samplePoints);
#else
      cell_ids          = initArray<vtkIntArray>("SegmentationId", samplePoints);
      cell_upNode       = initArray<vtkIntArray>("upNodeId", samplePoints);
      cell_downNode     = initArray<vtkIntArray>("downNodeId", samplePoints);
#endif
      point_regularMask = initArray<vtkCharArray>("RegularMask", samplePoints);
      point_scalar      = initArray<vtkFloatArray>("Scalar", samplePoints);

      if (params.advStats) {
         if (params.segm) {
#ifdef TTK_USE_64BIT_IDS
            cell_sizeArcs = initArray<vtkIdTypeArray>("RegionSize", samplePoints);
#else
            cell_sizeArcs = initArray<vtkIntArray>("RegionSize", samplePoints);
#endif
         }
         cell_spanArcs = initArray<vtkDoubleArray>("RegionSpan", samplePoints);
      }

      return 0;
   }

   inline bool hasPoint(const ttk::ftm::idVertex vertId)
   {
      return point_ids[vertId] != ttk::ftm::nullVertex;
   }

   inline void addPoint(const ttk::ftm::idVertex globalId, const ttk::SimplexId id, const float scalar,
                        const bool reg)
   {
      point_ids[globalId] = id;
      setPoint(id, scalar, reg);
   }

   inline void setPoint(const ttk::SimplexId id, const float scalar, const bool reg)
   {
      point_scalar->SetTuple1(id, scalar);
      point_regularMask->SetTuple1(id, reg);
   }

   inline void fillArrayCell(const ttk::SimplexId pos, const ttk::ftm::idSuperArc arcId,
                             LocalFTM& ftmTree, ttk::Triangulation* triangulation,
                             ttk::ftm::Params params)
   {
      const ttk::ftm::idNode idOffset = ftmTree.offset;
      ttk::ftm::FTMTree_MT*  tree     = ftmTree.tree.getTree(params.treeType);
      ttk::ftm::SuperArc*    arc      = tree->getSuperArc(arcId);

      if (params.normalize) {
         cell_ids->SetTuple1(pos, idOffset + arc->getNormalizedId());
      } else {
         cell_ids->SetTuple1(pos, idOffset + arcId);
      }

      cell_upNode->SetTuple1(pos, arc->getUpNodeId());
      cell_downNode->SetTuple1(pos, arc->getDownNodeId());

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs->SetTuple1(pos, tree->getArcSize(arcId));
         }

         float               downPoints[3];
         const ttk::ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ttk::ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ttk::ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ttk::ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         cell_spanArcs->SetTuple1(pos, ttk::Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkUnstructuredGrid* skeletonArcs, ttk::ftm::Params params)
   {
      // Some arcs might have been less sampled than the desired value, if they
      // have not enought regular vertices. Here we ensur that we will no keep
      // noise in these arrays.
      const size_t nbPoints = skeletonArcs->GetNumberOfPoints();
      const size_t nbCells  = skeletonArcs->GetNumberOfCells();

      cell_ids->SetNumberOfTuples(nbCells);
      skeletonArcs->GetCellData()->SetScalars(cell_ids);

      cell_upNode->SetNumberOfTuples(nbCells);
      cell_downNode->SetNumberOfTuples(nbCells);
      skeletonArcs->GetCellData()->AddArray(cell_upNode);
      skeletonArcs->GetCellData()->AddArray(cell_downNode);

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs->SetNumberOfTuples(nbCells);
            skeletonArcs->GetCellData()->AddArray(cell_sizeArcs);
         }
         cell_spanArcs->SetNumberOfTuples(nbCells);
         skeletonArcs->GetCellData()->AddArray(cell_spanArcs);
      }

      point_scalar->SetNumberOfTuples(nbPoints);
      skeletonArcs->GetPointData()->AddArray(point_scalar);
      point_regularMask->SetNumberOfTuples(nbPoints);
      skeletonArcs->GetPointData()->AddArray(point_regularMask);

      point_ids.clear();
   }
};

struct NodeData : public WrapperData{
#ifdef TTK_USE_64BIT_IDS
   vtkSmartPointer<vtkIdTypeArray>   ids;
   vtkSmartPointer<vtkIdTypeArray>   vertIds;
   vtkSmartPointer<vtkIdTypeArray>   regionSize;
   vtkSmartPointer<vtkIdTypeArray>   regionSpan;
#else
   vtkSmartPointer<vtkIntArray>   ids;
   vtkSmartPointer<vtkIntArray>   vertIds;
   vtkSmartPointer<vtkIntArray>   regionSize;
   vtkSmartPointer<vtkIntArray>   regionSpan;
#endif
   vtkSmartPointer<vtkFloatArray> scalars;
   vtkSmartPointer<vtkIntArray>   type;
   int scalarType;

   inline int init(std::vector<LocalFTM>& ftmTree, ttk::ftm::Params params)
   {
      ttk::ftm::idNode numberOfNodes = 0;
      for (auto& t : ftmTree) {
         ttk::ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfNodes += tree->getNumberOfNodes();
      }

#ifdef TTK_USE_64BIT_IDS
      ids     = initArray<vtkIdTypeArray>("NodeId", numberOfNodes);
      vertIds = initArray<vtkIdTypeArray>("VertexId", numberOfNodes);
#else
      ids     = initArray<vtkIntArray>("NodeId", numberOfNodes);
      vertIds = initArray<vtkIntArray>("VertexId", numberOfNodes);
#endif
      type    = initArray<vtkIntArray>("NodeType", numberOfNodes);
      scalars = initArray<vtkFloatArray>("Scalar", numberOfNodes);

      if (params.advStats) {
         if (params.segm) {
#ifdef TTK_USE_64BIT_IDS
            regionSize = initArray<vtkIdTypeArray>("RegionSize", numberOfNodes);
#else
            regionSize = initArray<vtkIntArray>("RegionSize", numberOfNodes);
#endif
         }
#ifdef TTK_USE_64BIT_IDS
         regionSpan = initArray<vtkIdTypeArray>("RegionSpan", numberOfNodes);
#else
         regionSpan = initArray<vtkIntArray>("RegionSpan", numberOfNodes);
#endif
      }

      return 0;
   }

   void setScalarType(const int s)
   {
      scalarType = s;
   }

   inline void fillArrayPoint(ttk::SimplexId arrIdx, const ttk::ftm::idNode nodeId, LocalFTM& ftmTree,
                              vtkDataArray* idMapper, ttk::Triangulation* triangulation,
                              ttk::ftm::Params params)
   {
      const ttk::ftm::idNode   idOffset   = ftmTree.offset;
      ttk::ftm::FTMTree_MT*    tree       = ftmTree.tree.getTree(params.treeType);
      const ttk::ftm::Node*    node       = tree->getNode(nodeId);
      // local (per cc) id
      const ttk::ftm::idVertex l_vertexId = node->getVertexId();
      // global id
      const ttk::ftm::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);
      float cellScalar = 0;
      switch (scalarType) {
         vtkTemplateMacro({ cellScalar = (float)tree->getValue<VTK_TT>(l_vertexId); });
      }

      ids->SetTuple1(arrIdx, idOffset + nodeId);
      scalars->SetTuple1(arrIdx, cellScalar);
      vertIds->SetTuple1(arrIdx, g_vertexId);
      type->SetTuple1(arrIdx, static_cast<int>(getNodeType(*tree, nodeId, params)));

      if (params.advStats) {
         ttk::ftm::idSuperArc saId = getAdjSa(node);
         if (params.segm) {
            regionSize->SetTuple1(arrIdx, tree->getArcSize(saId));
         }

         ttk::ftm::SuperArc* arc = tree->getSuperArc(saId);

         float               downPoints[3];
         const ttk::ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ttk::ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ttk::ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ttk::ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         regionSpan->SetTuple1(arrIdx, ttk::Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkPointData* pointData, ttk::ftm::Params params)
   {
      pointData->AddArray(ids);
      pointData->AddArray(scalars);
      pointData->AddArray(vertIds);
      pointData->SetScalars(type);
      if (params.advStats) {
         if (params.segm) {
            pointData->AddArray(regionSize);
         }
         pointData->AddArray(regionSpan);
      }
   }

  private:

   ttk::ftm::idSuperArc getAdjSa(const ttk::ftm::Node* node)
   {
      if (node->getNumberOfDownSuperArcs() == 1) {
         return node->getDownSuperArcId(0);
      }

      if (node->getNumberOfUpSuperArcs() == 1) {
         return node->getUpSuperArcId(0);
      }

      // Degenerate case, arbitrary choice
      if (node->getNumberOfDownSuperArcs()) {
         return node->getDownSuperArcId(0);
      }

      if (node->getNumberOfDownSuperArcs()) {
         return node->getDownSuperArcId(0);
      }

      // Empty node
#ifndef TTK_ENABLE_KAMIKAZE
      std::cerr << "[ttkFTMTree]: node without arcs:" << node->getVertexId() << std::endl;
#endif
      return ttk::ftm::nullSuperArc;
   }
};

struct VertData: public WrapperData {
#ifdef TTK_USE_64BIT_IDS
   vtkSmartPointer<vtkIdTypeArray>    ids;
   vtkSmartPointer<vtkIdTypeArray>    sizeRegion;
#else
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    sizeRegion;
#endif
   vtkSmartPointer<vtkDoubleArray> spanRegion;
   vtkSmartPointer<vtkCharArray>   typeRegion;

   inline int init(std::vector<LocalFTM>& ftmTrees, ttk::ftm::Params params)
   {
      if (!params.segm)
         return 0;

      ttk::ftm::idVertex numberOfVertices = 0;

      for (auto& t : ftmTrees) {
         ttk::ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfVertices += tree->getNumberOfVertices();
      }

#ifdef TTK_USE_64BIT_IDS
      ids        = initArray<vtkIdTypeArray>("SegmentationId", numberOfVertices);
#else
      ids        = initArray<vtkIntArray>("SegmentationId", numberOfVertices);
#endif
      typeRegion = initArray<vtkCharArray>("RegionType", numberOfVertices);

      if (params.advStats) {
#ifdef TTK_USE_64BIT_IDS
         sizeRegion = initArray<vtkIdTypeArray>("RegionSize", numberOfVertices);
#else
         sizeRegion = initArray<vtkIntArray>("RegionSize", numberOfVertices);
#endif
         spanRegion = initArray<vtkDoubleArray>("RegionSpan", numberOfVertices);
      }

      return 0;
   }

   void fillArrayPoint(const ttk::ftm::idSuperArc arcId, LocalFTM& l_tree,
                       ttk::Triangulation* triangulation, vtkDataArray* idMapper,
                       ttk::ftm::Params params)
   {
      if (!params.segm)
         return;

      ttk::ftm::FTMTree_MT* tree = l_tree.tree.getTree(params.treeType);
      const ttk::ftm::idNode idOffset = l_tree.offset;
      ttk::ftm::SuperArc* arc = tree->getSuperArc(arcId);

      const ttk::ftm::idNode   upNodeId     = arc->getUpNodeId();
      const ttk::ftm::Node*    upNode       = tree->getNode(upNodeId);
      const ttk::ftm::idVertex l_upVertexId = upNode->getVertexId();
      const ttk::ftm::idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
      const ttk::ftm::NodeType upNodeType   = getNodeType(*tree, upNodeId, params);
      float               coordUp[3];
      triangulation->getVertexPoint(l_upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const ttk::ftm::idNode   downNodeId     = arc->getDownNodeId();
      const ttk::ftm::Node*    downNode       = tree->getNode(downNodeId);
      const ttk::ftm::idVertex           l_downVertexId = downNode->getVertexId();
      const ttk::ftm::idVertex           g_downVertexId = idMapper->GetTuple1(l_downVertexId);
      const ttk::ftm::NodeType downNodeType   = getNodeType(*tree, downNodeId, params);
      float               coordDown[3];
      triangulation->getVertexPoint(l_downVertexId, coordDown[0], coordDown[1], coordDown[2]);

      const ttk::ftm::idVertex    regionSize = tree->getSuperArc(arcId)->getNumberOfRegularNodes();
      const double regionSpan = ttk::Geometry::distance(coordUp, coordDown);

      ttk::ftm::idSuperArc nid = arc->getNormalizedId();

      ttk::ftm::ArcType regionType;
      // RegionType
      if (upNodeType == ttk::ftm::NodeType::Local_minimum &&
          downNodeType == ttk::ftm::NodeType::Local_maximum)
      {
         regionType = ttk::ftm::ArcType::Min_arc;
      }
      else if (upNodeType == ttk::ftm::NodeType::Local_minimum ||
               downNodeType == ttk::ftm::NodeType::Local_minimum)
      {
         regionType = ttk::ftm::ArcType::Min_arc;
      }
      else if (upNodeType == ttk::ftm::NodeType::Local_maximum ||
               downNodeType == ttk::ftm::NodeType::Local_maximum)
      {
         regionType = ttk::ftm::ArcType::Max_arc;
      }
      else if (upNodeType == ttk::ftm::NodeType::Saddle1 &&
               downNodeType == ttk::ftm::NodeType::Saddle1)
      {
         regionType = ttk::ftm::ArcType::Saddle1_arc;
      }
      else if (upNodeType == ttk::ftm::NodeType::Saddle2 &&
               downNodeType == ttk::ftm::NodeType::Saddle2)
      {
         regionType = ttk::ftm::ArcType::Saddle2_arc;
      }
      else
      {
         regionType = ttk::ftm::ArcType::Saddle1_saddle2_arc;
      }

      // fill extrema and regular verts of this arc

      // critical points
      if(params.normalize){
         ids->SetTuple1(g_upVertexId   , idOffset + nid);
         ids->SetTuple1(g_downVertexId , idOffset + nid);
      } else{
         ids->SetTuple1(g_upVertexId   , idOffset + arcId);
         ids->SetTuple1(g_downVertexId , idOffset + arcId);
      }

      if (params.advStats) {
         sizeRegion->SetTuple1(g_upVertexId   , regionSize);
         sizeRegion->SetTuple1(g_downVertexId , regionSize);
         spanRegion->SetTuple1(g_upVertexId   , regionSpan);
         spanRegion->SetTuple1(g_downVertexId , regionSpan);
      }
      typeRegion->SetTuple1(g_upVertexId   , static_cast<char>(regionType));
      typeRegion->SetTuple1(g_downVertexId , static_cast<char>(regionType));

      // regular nodes
      for (const ttk::ftm::idVertex l_vertexId : *arc) {
         const ttk::ftm::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);
         if (params.normalize) {
            ids->SetTuple1(g_vertexId, idOffset + nid);
         } else {
            ids->SetTuple1(g_vertexId, idOffset + arcId);
         }
         if (params.advStats) {
            sizeRegion->SetTuple1(g_vertexId, regionSize);
            spanRegion->SetTuple1(g_vertexId, regionSpan);
         }
         typeRegion->SetTuple1(g_vertexId, static_cast<char>(regionType));
      }

   }

   void addArray(vtkPointData* pointData, ttk::ftm::Params params)
   {
      if (!params.segm)
         return;

      pointData->SetScalars(ids);

      if (params.advStats) {
         pointData->AddArray(sizeRegion);
         pointData->AddArray(spanRegion);
      }
      pointData->AddArray(typeRegion);
   }
};

#endif /* end of include guard: TTKFTMSTRUCTURES_H */
