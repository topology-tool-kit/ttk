/// \sa ttk::ftm::FTMTree
#ifndef _VTK_CONTOURFORESTS_H
#define _VTK_CONTOURFORESTS_H

// ttk code includes
#include <FTMTree.h>
#include <ttkWrapper.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

struct WrapperData {
   template<typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> initArray(const char* fieldName, size_t nbElmnt)
   {
      vtkSmartPointer<vtkArrayType> arr = vtkSmartPointer<vtkArrayType>::New();
      arr->SetName(fieldName);
      arr->SetNumberOfComponents(1);
      arr->SetNumberOfTuples(nbElmnt);

#ifndef TTK_ENABLE_KAMIKAZE
      if (!arr) {
         cerr << "[ttkFTMTree] Error, unable to allocate " << fieldName
              << " the program will likely crash" << endl;
      }
#endif
      return arr;
   }

   inline static ftm::NodeType getNodeType(ftm::FTMTree_MT* tree, const ftm::idNode nodeId,
                                           ftm::Params params)
   {
      const ftm::Node* node = tree->getNode(nodeId);
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
   vector<vtkIdType>               pointIds;
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    sizeArcs;
   vtkSmartPointer<vtkDoubleArray> spanArcs;
   vtkSmartPointer<vtkCharArray>   regularMask;

   inline int init(ftm::FTMTree_MT* tree, ftm::Params params)
   {
      pointIds.resize(tree->getNumberOfNodes(), ftm::nullNodes);

      const ftm::idSuperArc nbArcs  = tree->getNumberOfSuperArcs();
      const ftm::idSuperArc nbNodes = tree->getNumberOfNodes();
      const ftm::idSuperArc samplePoints = params.samplingLvl >= 0
                                               ? nbNodes + (nbArcs * params.samplingLvl)
                                               : tree->getNumberOfVertices();

      ids         = initArray<vtkIntArray>("SegmentationId", samplePoints);
      regularMask = initArray<vtkCharArray>("RegularMask", samplePoints);

      if(params.advStats){
         if (params.segm) {
            sizeArcs = initArray<vtkIntArray>("RegionSize", samplePoints);
         }
         spanArcs = initArray<vtkDoubleArray>("RegionSpan", samplePoints);
      }

      return 0;
   }

   inline void fillArray(const vtkIdType pos, const ftm::idSuperArc arcId, ftm::FTMTree_MT* tree,
                         Triangulation* triangulation, ftm::Params params)
   {
      ftm::SuperArc* arc = tree->getSuperArc(arcId);

      if (params.normalize) {
         ids->SetTuple1(pos, arc->getNormalizedId());
      } else {
         ids->SetTuple1(pos, arcId);
      }

      if (params.advStats) {
         if (params.segm) {
            sizeArcs->SetTuple1(pos, tree->getArcSize(arcId));
         }

         float               downPoints[3];
         const ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         spanArcs->SetTuple1(pos, Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkUnstructuredGrid* skeletonArcs, ftm::Params params)
   {
      // Some arcs might have been less sampled than the desired value, if they have not enought
      // regular vertices. Here we ensur that we will no keep noise in these arrays.
      const size_t nbPoints = skeletonArcs->GetNumberOfPoints();
      const size_t nbCells = nbPoints -1;

      ids->SetNumberOfTuples(nbCells);
      skeletonArcs->GetCellData()->AddArray(ids);

      if (params.advStats) {
         if (params.segm) {
            sizeArcs->SetNumberOfTuples(nbCells);
            skeletonArcs->GetCellData()->AddArray(sizeArcs);
         }
         spanArcs->SetNumberOfTuples(nbCells);
         skeletonArcs->GetCellData()->AddArray(spanArcs);
      }

      regularMask->SetNumberOfTuples(nbPoints);
      skeletonArcs->GetPointData()->AddArray(regularMask);

      pointIds.clear();
   }
};

struct NodeData : public WrapperData{
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkIntArray> vertIds;
   vtkSmartPointer<vtkIntArray> type;
   vtkSmartPointer<vtkIntArray> regionSize;
   vtkSmartPointer<vtkIntArray> regionSpan;

   inline int init(ftm::FTMTree_MT* tree, ftm::Params params)
   {
      const ftm::idNode numberOfNodes = tree->getNumberOfNodes();
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

   inline void fillArray(const ftm::idNode nodeId, ftm::FTMTree_MT* tree,
                         Triangulation* triangulation, ftm::Params params)
   {
      const ftm::Node*    node     = tree->getNode(nodeId);
      const ftm::idVertex vertexId = node->getVertexId();

      ids->SetTuple1(nodeId, nodeId);
      vertIds->SetTuple1(nodeId, vertexId);
      type->SetTuple1(nodeId, static_cast<int>(getNodeType(tree, nodeId, params)));

      if (params.advStats) {
         ftm::idSuperArc saId = getAdjSa(node);
         if (params.segm) {
            regionSize->SetTuple1(nodeId, tree->getArcSize(saId));
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

         regionSpan->SetTuple1(nodeId, Geometry::distance(downPoints, upPoints));
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

   int init(ftm::FTMTree_MT* tree, ftm::Params params)
   {
      if (!params.segm)
         return 0;

      const ftm::idVertex numberOfVertices = tree->getNumberOfVertices();

      ids        = initArray<vtkIntArray>("SegmentationId", numberOfVertices);
      typeRegion = initArray<vtkCharArray>("RegionType", numberOfVertices);

      if (params.advStats){
         sizeRegion = initArray<vtkIntArray>("RegionSize", numberOfVertices);
         spanRegion = initArray<vtkDoubleArray>("RegionSpan", numberOfVertices);
      }

      return 0;
   }

   void fillArray(const ftm::idSuperArc arcId, ftm::FTMTree_MT* tree, Triangulation* triangulation,
                  ftm::Params params)
   {
      if (!params.segm)
         return;

      ftm::SuperArc* arc = tree->getSuperArc(arcId);

      const int           upNodeId   = arc->getUpNodeId();
      const ftm::Node*    upNode     = tree->getNode(upNodeId);
      const int           upVertexId = upNode->getVertexId();
      const ftm::NodeType upNodeType = getNodeType(tree, upNodeId, params);
      float               coordUp[3];
      triangulation->getVertexPoint(upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const int           downNodeId   = arc->getDownNodeId();
      const ftm::Node*    downNode     = tree->getNode(downNodeId);
      const int           downVertexId = downNode->getVertexId();
      const ftm::NodeType downNodeType = getNodeType(tree, downNodeId, params);
      float               coordDown[3];
      triangulation->getVertexPoint(downVertexId, coordDown[0], coordDown[1], coordDown[2]);

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
         ids->SetTuple1(upVertexId   , nid);
         ids->SetTuple1(downVertexId , nid);
      } else{
         ids->SetTuple1(upVertexId   , arcId);
         ids->SetTuple1(downVertexId , arcId);
      }

      if (params.advStats) {
         sizeRegion->SetTuple1(upVertexId   , regionSize);
         sizeRegion->SetTuple1(downVertexId , regionSize);
         spanRegion->SetTuple1(upVertexId   , regionSpan);
         spanRegion->SetTuple1(downVertexId , regionSpan);
      }
      typeRegion->SetTuple1(upVertexId   , static_cast<char>(regionType));
      typeRegion->SetTuple1(downVertexId , static_cast<char>(regionType));

      // regular nodes
      for (const ftm::idVertex vertexId : *arc) {
         if (params.normalize) {
            ids->SetTuple1(vertexId, nid);
         } else {
            ids->SetTuple1(vertexId, arcId);
         }
         if (params.advStats) {
            sizeRegion->SetTuple1(vertexId, regionSize);
            spanRegion->SetTuple1(vertexId, regionSpan);
         }
         typeRegion->SetTuple1(vertexId, static_cast<char>(regionType));
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


class VTKFILTERSCORE_EXPORT ttkFTMTree : public vtkDataSetAlgorithm, public Wrapper
{
  public:
   static ttkFTMTree* New();

   vtkTypeMacro(ttkFTMTree, vtkDataSetAlgorithm);

   // default ttk setters
   vtkSetMacro(debugLevel_, int);

   void SetThreadNumber(int threadNumber)
   {
      ThreadNumber = threadNumber;
      SetThreads();
   }

   void SetUseAllCores(bool onOff)
   {
      UseAllCores = onOff;
      SetThreads();
   }
   // end of default ttk setters

   vtkSetMacro(ScalarField, string);
   vtkGetMacro(ScalarField, string);

   vtkSetMacro(UseInputOffsetScalarField, int);
   vtkGetMacro(UseInputOffsetScalarField, int);

   vtkSetMacro(InputOffsetScalarFieldName, string);
   vtkGetMacro(InputOffsetScalarFieldName, string);

   vtkSetMacro(ScalarFieldId, int);
   vtkGetMacro(ScalarFieldId, int);

   vtkSetMacro(OffsetFieldId, int);
   vtkGetMacro(OffsetFieldId, int);

   // Parameters uses a structure, we can't use vtkMacro on them
   void SetTreeType(const int type)
   {
       params_.treeType = (ftm::TreeType)type;
       Modified();
   }

   ftm::TreeType GetTreeType(void) const
   {
       return params_.treeType;
   }

   void SetWithSegmentation(const bool segm)
   {
      params_.segm = segm;
      Modified();
   }

   bool GetWithSegmentation(void) const
   {
      return params_.segm;
   }

   void SetWithNormalize(const bool norm)
   {
      params_.normalize = norm;
      Modified();
   }

   bool GetWithNormalize(void) const
   {
      return params_.normalize;
   }

   void SetWithAdvStats(const bool adv)
   {
      params_.advStats = adv;
      Modified();
   }

   bool GetWithAdvStats(void) const
   {
      return params_.advStats;
   }

   void SetSuperArcSamplingLevel(int lvl)
   {
      params_.samplingLvl = lvl;
      Modified();
   }

   int GetSuperArcSamplingLevel(void) const
   {
      return params_.samplingLvl;
   }

   int setupTriangulation(vtkDataSet* input);
   int getScalars(vtkDataSet* input);
   int getOffsets(vtkDataSet* input);

   int getSkeletonNodes(ftm::FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(ftm::FTMTree_MT* tree,
                            ftm::idSuperArc arcId,
                            vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs,
                            ArcData& arcData);

   int addSampledSkeletonArc(ftm::FTMTree_MT* tree,
                             ftm::idSuperArc arcId,
                             const int samplingLevel,
                             vtkPoints* points,
                             vtkUnstructuredGrid* skeletonArcs,
                             ArcData& arcData);

   int addCompleteSkeletonArc(ftm::FTMTree_MT* tree,
                              ftm::idSuperArc arcId,
                              vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs,
                              ArcData& arcData);

   int getSkeletonArcs(ftm::FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonArcs);

   int getSegmentation(ftm::FTMTree_MT* tree, vtkDataSet* input, vtkDataSet* outputSegmentation);

#ifdef withStatsTime
   void printCSVStats();
   void printCSVTree(const ftm::FTMTree_MT* const tree) const;
#endif

  protected:
   ttkFTMTree();
   ~ttkFTMTree();

   TTK_SETUP();

   virtual int FillInputPortInformation(int port, vtkInformation* info);
   virtual int FillOutputPortInformation(int port, vtkInformation* info);

  private:
   string ScalarField;
   bool   UseInputOffsetScalarField;
   string InputOffsetScalarFieldName;
   int    ScalarFieldId;
   int    OffsetFieldId;

   ftm::Params params_;

   Triangulation* triangulation_;
   ftm::FTMTree   ftmTree_;
   vtkDataArray*  inputScalars_;
   vtkIntArray*   offsets_;
   vtkDataArray*  inputOffsets_;
   bool           hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
