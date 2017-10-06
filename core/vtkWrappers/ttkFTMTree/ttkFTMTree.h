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

struct ArcData {
   vector<vtkIdType>               pointIds;
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    normalizedIds;
   vtkSmartPointer<vtkIntArray>    sizeArcs;
   vtkSmartPointer<vtkDoubleArray> spanArcs;
   vtkSmartPointer<vtkCharArray>   regularMask;
#ifdef withStatsTime
   vtkSmartPointer<vtkFloatArray> startArcs;
   vtkSmartPointer<vtkFloatArray> endArcs;
   vtkSmartPointer<vtkFloatArray> timeArcs;
   vtkSmartPointer<vtkIntArray>   origArcs;
   vtkSmartPointer<vtkIntArray>   tasksArcs;
#endif

   int init(const ftm::idNode nbNodes, ftm::Params params)
   {
      pointIds.resize(nbNodes, ftm::nullNodes);

      if(params.normalize){
         normalizedIds = vtkSmartPointer<vtkIntArray>::New();
         normalizedIds->SetName("SegmentationId");
         normalizedIds->SetNumberOfComponents(1);
      } else {
         ids = vtkSmartPointer<vtkIntArray>::New();
         ids->SetName("SegmentationId");
         ids->SetNumberOfComponents(1);
      }

      sizeArcs = vtkSmartPointer<vtkIntArray>::New();
      sizeArcs->SetName("RegionSize");
      sizeArcs->SetNumberOfComponents(1);

      spanArcs = vtkSmartPointer<vtkDoubleArray>::New();
      spanArcs->SetName("RegionSpan");
      spanArcs->SetNumberOfComponents(1);

      regularMask = vtkSmartPointer<vtkCharArray>::New();
      regularMask->SetName("RegularMask");
      regularMask->SetNumberOfComponents(1);

#ifdef withStatsTime
      startArcs = vtkSmartPointer<vtkFloatArray>::New();
      startArcs->SetName("Start");
      startArcs->SetNumberOfComponents(1);

      endArcs = vtkSmartPointer<vtkFloatArray>::New();
      endArcs->SetName("End");
      endArcs->SetNumberOfComponents(1);

      timeArcs = vtkSmartPointer<vtkFloatArray>::New();
      timeArcs->SetName("Time");
      timeArcs->SetNumberOfComponents(1);

      origArcs = vtkSmartPointer<vtkIntArray>::New();
      origArcs->SetName("Origin");
      origArcs->SetNumberOfComponents(1);

      tasksArcs = vtkSmartPointer<vtkIntArray>::New();
      tasksArcs->SetName("Tasks");
      tasksArcs->SetNumberOfComponents(1);
#endif

// Check
#ifndef withKamikaze

      if (!ids && !normalizedIds) {
         cerr << "[ttkFTMTree] Error : vtkIntArray ids arcs allocation problem." << endl;
         return -1;
      }

      if (!sizeArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray size allocation problem." << endl;
         return -1;
      }

      if (!spanArcs) {
         cerr << "[ttkFTMTree] Error : vtkDoubleArray span allocation problem." << endl;
         return -1;
      }

      if (!regularMask) {
         cerr << "[ttkFTMTree] Error : vtkCharArray span allocation problem." << endl;
         return -1;
      }

# ifdef withStatsTime

      if (!startArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray start allocation problem." << endl;
         return -1;
      }

      if (!endArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray end allocation problem." << endl;
         return -1;
      }

      if (!timeArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray time allocation problem." << endl;
         return -1;
      }

      if (!origArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray origin allocation problem." << endl;
         return -1;
      }

      if (!tasksArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray tasks allocation problem." << endl;
         return -1;
      }

# endif
#endif
      return 0;
   }

   void fillArray(const ftm::idSuperArc arcId, Triangulation* triangulation, ftm::FTMTree_MT* tree,
                  ftm::Params params, bool reg=false)
   {
      ftm::SuperArc* arc = tree->getSuperArc(arcId);

      float          downPoints[3];
      const ftm::idVertex downNodeId   = tree->getLowerNodeId(arc);
      const ftm::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
      triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

      float          upPoints[3];
      const ftm::idVertex upNodeId   = tree->getUpperNodeId(arc);
      const ftm::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
      triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

      if (params.normalize) {
         normalizedIds->InsertNextTuple1(arc->getNormalizedId());
      } else {
         ids->InsertNextTuple1(arcId);
      }

      if (params.advStats) {
         if (params.segm) {
            sizeArcs->InsertNextTuple1(tree->getArcSize(arcId));
         }
         spanArcs->InsertNextTuple1(Geometry::distance(downPoints, upPoints));
      }

#ifdef withStatsTime
      startArcs->InsertNextTuple1(tree->getArcStart(arcId));
      endArcs  ->InsertNextTuple1(tree->getArcEnd(arcId));
      timeArcs ->InsertNextTuple1(tree->getArcEnd(arcId) - tree->getArcStart(arcId));
      origArcs ->InsertNextTuple1(tree->getArcOrig(arcId));
      tasksArcs->InsertNextTuple1(tree->getArcActiveTasks(arcId));
#endif
   }

   void addArray(vtkUnstructuredGrid *skeletonArcs, ftm::Params params){
      if (params.normalize) {
         skeletonArcs->GetCellData()->AddArray(normalizedIds);
      } else {
         skeletonArcs->GetCellData()->AddArray(ids);
      }

      if (params.advStats) {
         if (params.segm) {
            skeletonArcs->GetCellData()->AddArray(sizeArcs);
         }
         skeletonArcs->GetCellData()->AddArray(spanArcs);
      }

      skeletonArcs->GetPointData()->AddArray(regularMask);
#ifdef withStatsTime
      skeletonArcs->GetCellData()->AddArray(startArcs);
      skeletonArcs->GetCellData()->AddArray(endArcs);
      skeletonArcs->GetCellData()->AddArray(timeArcs);
      skeletonArcs->GetCellData()->AddArray(origArcs);
      skeletonArcs->GetCellData()->AddArray(tasksArcs);
#endif
      pointIds.clear();
   }
};

struct NodeData {
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkIntArray>   vertIds;
   vtkSmartPointer<vtkIntArray> type;
   vtkSmartPointer<vtkIntArray> regionSize;
   vtkSmartPointer<vtkIntArray> regionSpan;

   int init(const ftm::idNode numberOfNodes, ftm::Params params)
   {
      ids = vtkSmartPointer<vtkIntArray>::New();
      ids->SetName("NodeId");
      ids->SetNumberOfComponents(1);

      vertIds = vtkSmartPointer<vtkIntArray>::New();
      vertIds->SetName("VertexId");
      vertIds->SetNumberOfComponents(1);

      type = vtkSmartPointer<vtkIntArray>::New();
      type->SetName("NodeType");
      type->SetNumberOfComponents(1);

      regionSize = vtkSmartPointer<vtkIntArray>::New();
      regionSize->SetName("RegionSize");
      regionSize->SetNumberOfComponents(1);

      regionSpan = vtkSmartPointer<vtkIntArray>::New();
      regionSpan->SetName("RegionSpan");
      regionSpan->SetNumberOfComponents(1);
// Check
#ifndef withKamikaze

      if (!ids) {
         cerr << "[ttkFTMTree] Error : vtkIntArray ids nodes allocation problem." << endl;
         return -1;
      }

      if (!vertIds) {
         cerr << "[ttkFTMTree] Error : vtkIntArray verts nodes allocation problem." << endl;
         return -1;
      }

      if (!type) {
         cerr << "[ttkFTMTree] Error : vtkIntArray type nodes allocation problem." << endl;
         return -1;
      }

      if (!regionSize) {
         cerr << "[ttkFTMTree] Error : vtkIntArray regionSize nodes allocation problem." << endl;
         return -1;
      }

      if (!regionSpan) {
         cerr << "[ttkFTMTree] Error : vtkIntArray regionSpan nodes allocation problem." << endl;
         return -1;
      }
#endif

   }

   void fillArray(const ftm::idNode nodeId, Triangulation* triangulation, ftm::FTMTree_MT* tree,
                  ftm::Params params)
   {
      const ftm::Node* node = tree->getNode(nodeId);
      const ftm::idVertex vertexId = node->getVertexId();
      ids->InsertNextTuple1(nodeId);
      vertIds->InsertNextTuple1(vertexId);
      type->InsertNextTuple1(static_cast<int>(getNodeType(node, params)));
      regionSize->InsertNextTuple1(0);
      regionSpan->InsertNextTuple1(0);
      cout << "region span in nodes, TODO" << endl;
   }

   void addArray(vtkPointData* pointData, ftm::Params params)
   {
      pointData->AddArray(ids);
      pointData->AddArray(vertIds);
      pointData->AddArray(type);
      if (params.advStats) {
         pointData->AddArray(regionSize);
         pointData->AddArray(regionSpan);
      }
   }

   static ftm::NodeType getNodeType(const ftm::Node* node, ftm::Params params)
   {
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

//TODO Unify init and avoid alloc when uneeded
//TODO inline function + static

struct VertData {
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    normalizedIdVert;
   vtkSmartPointer<vtkIntArray>    sizeRegion;
   vtkSmartPointer<vtkDoubleArray> spanRegion;
   vtkSmartPointer<vtkCharArray>   typeRegion;

   int init(const ftm::idVertex numberOfVertices, ftm::Params params)
   {
      if (!params.segm)
         return 0;

      if (params.normalize) {
         normalizedIdVert = vtkSmartPointer<vtkIntArray>::New();
         normalizedIdVert->SetName("SegmentationId");
         normalizedIdVert->SetNumberOfComponents(1);
         normalizedIdVert->SetNumberOfTuples(numberOfVertices);
      } else {
         ids = vtkSmartPointer<vtkIntArray>::New();
         ids->SetName("SegmentationId");
         ids->SetNumberOfComponents(1);
         ids->SetNumberOfTuples(numberOfVertices);
      }

      if (params.advStats){
         sizeRegion = vtkSmartPointer<vtkIntArray>::New();
         sizeRegion->SetName("RegionSize");
         sizeRegion->SetNumberOfComponents(1);
         sizeRegion->SetNumberOfTuples(numberOfVertices);

         spanRegion = vtkSmartPointer<vtkDoubleArray>::New();
         spanRegion->SetName("RegionSpan");
         spanRegion->SetNumberOfComponents(1);
         spanRegion->SetNumberOfTuples(numberOfVertices);
      }

      typeRegion = vtkSmartPointer<vtkCharArray>::New();
      typeRegion->SetName("RegionType");
      typeRegion->SetNumberOfComponents(1);
      typeRegion->SetNumberOfTuples(numberOfVertices);

// Check
#ifndef withKamikaze

      if (!ids && !normalizedIdVert) {
         cerr << "[ttkFTMTree] Error : vtkIntArray ids verts allocation problem." << endl;
         return -2;
      }

      if (params.advStats) {
         if (!sizeRegion) {
            cerr << "[ttkFTMTree] Error : vtkIntArray region size allocation problem." << endl;
            return -2;
         }

         if (!spanRegion) {
            cerr << "[ttkFTMTree] Error : vtkDoubleArray region span allocation problem." << endl;
            return -2;
         }
      }

      if (!typeRegion) {
         cerr << "[ttkFTMTree] Error : vtkCharArray region type allocation problem." << endl;
         return -2;
      }

#endif
      return 0;
   }

   void fillArray(const ftm::idSuperArc arcId, Triangulation* triangulation, ftm::FTMTree_MT* tree,
                  ftm::Params params)
   {
      if (!params.segm)
         return;

      auto getNodeType = [&](const ftm::idNode nodeId) {
         ftm::Node* node = tree->getNode(nodeId);
         int upDegree{};
         int downDegree{};
         if (tree->isST()) {
            downDegree = node->getNumberOfUpSuperArcs();
            upDegree   = node->getNumberOfDownSuperArcs();
         } else {
            upDegree   = node->getNumberOfUpSuperArcs();
            downDegree = node->getNumberOfDownSuperArcs();
         }
         int degree = upDegree + downDegree;

         // saddle point
         if (degree > 1) {
            if (upDegree == 2 and downDegree == 1)
               return ftm::NodeType::Saddle2;
            else if (upDegree == 1 && downDegree == 2)
               return ftm::NodeType::Saddle1;
            else if (upDegree == 1 && downDegree == 1)
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
      };

      ftm::SuperArc* arc = tree->getSuperArc(arcId);

      const int           upNodeId   = arc->getUpNodeId();
      const ftm::Node*    upNode     = tree->getNode(upNodeId);
      const int           upVertexId = upNode->getVertexId();
      const ftm::NodeType upNodeType = getNodeType(upNodeId);
      float               coordUp[3];
      triangulation->getVertexPoint(upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const int           downNodeId   = arc->getDownNodeId();
      const ftm::Node*    downNode     = tree->getNode(downNodeId);
      const int           downVertexId = downNode->getVertexId();
      const ftm::NodeType downNodeType = getNodeType(downNodeId);
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
         normalizedIdVert->SetTuple1(upVertexId   , nid);
         normalizedIdVert->SetTuple1(downVertexId , nid);
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
            normalizedIdVert->SetTuple1(vertexId, nid);
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

   void addArray(vtkPointData* pointData, ftm::Params params){
      if (!params.segm)
         return;

      if(params.normalize){
         pointData->AddArray(normalizedIdVert);
      } else {
         pointData->AddArray(ids);
      }

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

   vtkSetMacro(SuperArcSamplingLevel, int);
   vtkGetMacro(SuperArcSamplingLevel, int);

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

   bool SetWithAdvStats(const bool adv)
   {
      params_.advStats = adv;
      Modified();
   }

   bool GetWithAdvStats(void) const
   {
      return params_.advStats;
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
   int    SuperArcSamplingLevel;

   ftm::Params params_;

   Triangulation* triangulation_;
   ftm::FTMTree        ftmTree_;
   vtkDataArray*  inputScalars_;
   vtkIntArray*   offsets_;
   vtkDataArray*  inputOffsets_;
   bool           hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
