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
   ftm::FTMTree tree;
   ftm::idNode  offset;
};

struct WrapperData {
   template<typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> initArray(const char* fieldName, size_t nbElmnt)
   {
      vtkSmartPointer<vtkArrayType> arr = vtkSmartPointer<vtkArrayType>::New();
      arr->SetName(fieldName);
      arr->SetNumberOfComponents(1);
      arr->SetNumberOfTuples(nbElmnt);

#ifndef withKamikaze
      if (!arr) {
         cerr << "[ttkFTMTree] Error, unable to allocate " << fieldName
              << " the program will likely crash" << endl;
      }
#endif
      return arr;
   }

   inline static ftm::NodeType getNodeType(ftm::FTMTree_MT& tree, const ftm::idNode nodeId,
                                           ftm::Params params)
   {
      const ftm::Node* node = tree.getNode(nodeId);
      int upDegree{};
      int downDegree{};
      if (params.treeType == ftm::TreeType::Join or params.treeType == ftm::TreeType::Contour) {
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
            return ftm::NodeType::Saddle2;
         else if (upDegree == 1 and downDegree == 2)
            return ftm::NodeType::Saddle1;
         else if (upDegree == 1 and downDegree == 1)
            return ftm::NodeType::Regular;
         else
            return ftm::NodeType::Degenerate;
      }
      // local extremum
      else {
         if (upDegree)
            return ftm::NodeType::Local_minimum;
         else
            return ftm::NodeType::Local_maximum;
      }
   }
};

struct ArcData : public WrapperData {
   vector<vtkIdType>               point_ids;
   vtkSmartPointer<vtkCharArray>   point_regularMask;
   vtkSmartPointer<vtkFloatArray>  point_scalar;
   vtkSmartPointer<vtkIntArray>    cell_ids;
   vtkSmartPointer<vtkIntArray>    cell_sizeArcs;
   vtkSmartPointer<vtkDoubleArray> cell_spanArcs;

   inline int init(vector<LocalFTM>& ftmTree, ftm::Params params)
   {
      ftm::idSuperArc nbArcs       = 0;
      ftm::idSuperArc nbNodes      = 0;
      ftm::idSuperArc samplePoints = 0;
      ftm::idVertex   nbVerts      = 0;

      for (auto& t : ftmTree) {
         ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         nbArcs += tree->getNumberOfSuperArcs();
         nbNodes += tree->getNumberOfNodes();
         samplePoints += params.samplingLvl >= 0
                             ? tree->getNumberOfNodes() + (nbArcs * params.samplingLvl)
                             : tree->getNumberOfVertices();
         nbVerts += tree->getNumberOfVertices();
      }

      point_ids.resize(nbVerts, ftm::nullVertex);
      cell_ids          = initArray<vtkIntArray>("SegmentationId", samplePoints);
      point_regularMask = initArray<vtkCharArray>("RegularMask", samplePoints);
      point_scalar      = initArray<vtkFloatArray>("Scalar", samplePoints);

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs = initArray<vtkIntArray>("RegionSize", samplePoints);
         }
         cell_spanArcs = initArray<vtkDoubleArray>("RegionSpan", samplePoints);
      }

      return 0;
   }

   inline bool hasPoint(const ftm::idVertex vertId)
   {
      return point_ids[vertId] != ftm::nullVertex;
   }

   inline void addPoint(const ftm::idVertex globalId, const vtkIdType id, const float scalar,
                        const bool reg)
   {
      point_ids[globalId] = id;
      setPoint(id, scalar, reg);
   }

   inline void setPoint(const vtkIdType id, const float scalar, const bool reg)
   {
      point_scalar->SetTuple1(id, scalar);
      point_regularMask->SetTuple1(id, reg);
   }

   inline void fillArrayCell(const vtkIdType pos, const ftm::idSuperArc arcId, LocalFTM& ftmTree,
                             Triangulation* triangulation, ftm::Params params)
   {
      const ftm::idNode idOffset = ftmTree.offset;
      ftm::FTMTree_MT*  tree     = ftmTree.tree.getTree(params.treeType);
      ftm::SuperArc*    arc      = tree->getSuperArc(arcId);

      if (params.normalize) {
         cell_ids->SetTuple1(pos, idOffset + arc->getNormalizedId());
      } else {
         cell_ids->SetTuple1(pos, idOffset + arcId);
      }

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs->SetTuple1(pos, tree->getArcSize(arcId));
         }

         float               downPoints[3];
         const ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         cell_spanArcs->SetTuple1(pos, Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkUnstructuredGrid* skeletonArcs, ftm::Params params)
   {
      // Some arcs might have been less sampled than the desired value, if they have not enought
      // regular vertices. Here we ensur that we will no keep noise in these arrays.
      const size_t nbPoints = skeletonArcs->GetNumberOfPoints();
      const size_t nbCells  = skeletonArcs->GetNumberOfCells();

      cell_ids->SetNumberOfTuples(nbCells);
      skeletonArcs->GetCellData()->AddArray(cell_ids);

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
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkIntArray> vertIds;
   vtkSmartPointer<vtkIntArray> type;
   vtkSmartPointer<vtkIntArray> regionSize;
   vtkSmartPointer<vtkIntArray> regionSpan;

   inline int init(vector<LocalFTM>& ftmTree, ftm::Params params)
   {
      ftm::idNode numberOfNodes = 0;
      for (auto& t : ftmTree) {
         ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfNodes += tree->getNumberOfNodes();
      }

      ids     = initArray<vtkIntArray>("NodeId", numberOfNodes);
      vertIds = initArray<vtkIntArray>("VertexId", numberOfNodes);
      type    = initArray<vtkIntArray>("NodeType", numberOfNodes);

      if (params.advStats) {
         if (params.segm) {
            regionSize = initArray<vtkIntArray>("RegionSize", numberOfNodes);
         }
         regionSpan = initArray<vtkIntArray>("RegionSpan", numberOfNodes);
      }

      return 0;
   }

   inline void fillArrayPoint(vtkIdType arrIdx, const ftm::idNode nodeId, LocalFTM& ftmTree,
                              vtkDataArray* idMapper, Triangulation* triangulation,
                              ftm::Params params)
   {
      const ftm::idNode   idOffset   = ftmTree.offset;
      ftm::FTMTree_MT*    tree       = ftmTree.tree.getTree(params.treeType);
      const ftm::Node*    node       = tree->getNode(nodeId);
      const ftm::idVertex l_vertexId = node->getVertexId();
      const ftm::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);

      ids->SetTuple1(arrIdx, idOffset + nodeId);
      vertIds->SetTuple1(arrIdx, g_vertexId);
      type->SetTuple1(arrIdx, static_cast<int>(getNodeType(*tree, nodeId, params)));

      if (params.advStats) {
         ftm::idSuperArc saId = getAdjSa(node);
         if (params.segm) {
            regionSize->SetTuple1(arrIdx, tree->getArcSize(saId));
         }

         ftm::SuperArc* arc = tree->getSuperArc(saId);

         float               downPoints[3];
         const ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         regionSpan->SetTuple1(arrIdx, Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkPointData* pointData, ftm::Params params)
   {
      pointData->AddArray(ids);
      pointData->AddArray(vertIds);
      pointData->AddArray(type);
      if (params.advStats) {
         if (params.segm) {
            pointData->AddArray(regionSize);
         }
         pointData->AddArray(regionSpan);
      }
   }

  private:

   ftm::idSuperArc getAdjSa(const ftm::Node* node)
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
      cerr << "[ttkFTMTree]: node without arcs:" << node->getVertexId() << endl;
      return ftm::nullSuperArc;
   }
};

struct VertData: public WrapperData {
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    sizeRegion;
   vtkSmartPointer<vtkDoubleArray> spanRegion;
   vtkSmartPointer<vtkCharArray>   typeRegion;

   inline int init(vector<LocalFTM>& ftmTrees, ftm::Params params)
   {
      if (!params.segm)
         return 0;

      ftm::idVertex numberOfVertices = 0;

      for (auto& t : ftmTrees) {
         ftm::FTMTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfVertices += tree->getNumberOfVertices();
      }

      ids        = initArray<vtkIntArray>("SegmentationId", numberOfVertices);
      typeRegion = initArray<vtkCharArray>("RegionType", numberOfVertices);

      if (params.advStats) {
         sizeRegion = initArray<vtkIntArray>("RegionSize", numberOfVertices);
         spanRegion = initArray<vtkDoubleArray>("RegionSpan", numberOfVertices);
      }

      return 0;
   }

   void fillArrayPoint(const ftm::idSuperArc arcId, LocalFTM& l_tree, Triangulation* triangulation,
                       vtkDataArray* idMapper, ftm::Params params)
   {
      if (!params.segm)
         return;

      ftm::FTMTree_MT* tree = l_tree.tree.getTree(params.treeType);
      const ftm::idNode idOffset = l_tree.offset;
      ftm::SuperArc* arc = tree->getSuperArc(arcId);

      const ftm::idNode   upNodeId     = arc->getUpNodeId();
      const ftm::Node*    upNode       = tree->getNode(upNodeId);
      const ftm::idVertex l_upVertexId = upNode->getVertexId();
      const ftm::idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
      const ftm::NodeType upNodeType   = getNodeType(*tree, upNodeId, params);
      float               coordUp[3];
      triangulation->getVertexPoint(l_upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const ftm::idNode   downNodeId     = arc->getDownNodeId();
      const ftm::Node*    downNode       = tree->getNode(downNodeId);
      const int           l_downVertexId = downNode->getVertexId();
      const int           g_downVertexId = idMapper->GetTuple1(l_downVertexId);
      const ftm::NodeType downNodeType   = getNodeType(*tree, downNodeId, params);
      float               coordDown[3];
      triangulation->getVertexPoint(l_downVertexId, coordDown[0], coordDown[1], coordDown[2]);

      const int    regionSize = tree->getSuperArc(arcId)->getNumberOfRegularNodes();
      const double regionSpan = Geometry::distance(coordUp, coordDown);

      ftm::idSuperArc nid = arc->getNormalizedId();

      ftm::ArcType regionType;
      // RegionType
      if (upNodeType == ftm::NodeType::Local_minimum && downNodeType == ftm::NodeType::Local_maximum)
         regionType = ftm::ArcType::Min_arc;
      else if (upNodeType == ftm::NodeType::Local_minimum || downNodeType == ftm::NodeType::Local_minimum)
         regionType = ftm::ArcType::Min_arc;
      else if (upNodeType == ftm::NodeType::Local_maximum || downNodeType == ftm::NodeType::Local_maximum)
         regionType = ftm::ArcType::Max_arc;
      else if (upNodeType == ftm::NodeType::Saddle1 && downNodeType == ftm::NodeType::Saddle1)
         regionType = ftm::ArcType::Saddle1_arc;
      else if (upNodeType == ftm::NodeType::Saddle2 && downNodeType == ftm::NodeType::Saddle2)
         regionType = ftm::ArcType::Saddle2_arc;
      else
         regionType = ftm::ArcType::Saddle1_saddle2_arc;

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
      for (const ftm::idVertex l_vertexId : *arc) {
         const ftm::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);
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

   void addArray(vtkPointData* pointData, ftm::Params params)
   {
      if (!params.segm)
         return;

      pointData->AddArray(ids);

      if (params.advStats) {
         pointData->AddArray(sizeRegion);
         pointData->AddArray(spanRegion);
      }
      pointData->AddArray(typeRegion);
   }
};

#endif /* end of include guard: TTKFTMSTRUCTURES_H */
