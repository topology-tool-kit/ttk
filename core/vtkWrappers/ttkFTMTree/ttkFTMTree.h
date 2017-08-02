/// \sa ttk::FTMTree
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
   vtkSmartPointer<vtkIntArray>    idArcs;
   vtkSmartPointer<vtkIntArray>    normalizedIdArc;
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

   int init(const idNode nbNodes, Params params)
   {
      pointIds.resize(nbNodes, nullNodes);

      if(params.normalize){
         normalizedIdArc = vtkSmartPointer<vtkIntArray>::New();
         normalizedIdArc->SetName("SegmentationId");
         normalizedIdArc->SetNumberOfComponents(1);
      } else {
         idArcs = vtkSmartPointer<vtkIntArray>::New();
         idArcs->SetName("SegmentationId");
         idArcs->SetNumberOfComponents(1);
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

      if (!idArcs && !normalizedIdArc) {
         cerr << "[ttkFTMTree] Error : vtkIntArray id allocation problem." << endl;
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

   void fillArray(const idSuperArc arcId, Triangulation* triangulation, FTMTree_MT* tree,
                  Params params, bool reg=false)
   {
      SuperArc* arc = tree->getSuperArc(arcId);

      float          downPoints[3];
      const idVertex downNodeId   = tree->getLowerNodeId(arc);
      const idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
      triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

      float          upPoints[3];
      const idVertex upNodeId   = tree->getUpperNodeId(arc);
      const idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
      triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

      if (params.normalize) {
         normalizedIdArc->InsertNextTuple1(arc->getNormalizedId());
      } else {
         idArcs->InsertNextTuple1(arcId);
      }
      if (params.segm) {
         sizeArcs->InsertNextTuple1(tree->getArcSize(arcId));
      }
      spanArcs->InsertNextTuple1(Geometry::distance(downPoints, upPoints));
#ifdef withStatsTime
      startArcs->InsertNextTuple1(tree->getArcStart(arcId));
      endArcs  ->InsertNextTuple1(tree->getArcEnd(arcId));
      timeArcs ->InsertNextTuple1(tree->getArcEnd(arcId) - tree->getArcStart(arcId));
      origArcs ->InsertNextTuple1(tree->getArcOrig(arcId));
      tasksArcs->InsertNextTuple1(tree->getArcActiveTasks(arcId));
#endif
   }

   void addArray(vtkUnstructuredGrid *skeletonArcs, Params params){
      if (params.normalize) {
         skeletonArcs->GetCellData()->AddArray(normalizedIdArc);
      } else {
         skeletonArcs->GetCellData()->AddArray(idArcs);
      }
      if(params.segm)
      {
         skeletonArcs->GetCellData()->AddArray(sizeArcs);
      }
      skeletonArcs->GetCellData()->AddArray(spanArcs);
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

struct VertData {
   vtkSmartPointer<vtkIntArray>    idVerts;
   vtkSmartPointer<vtkIntArray>    normalizedIdVert;
   vtkSmartPointer<vtkIntArray>    sizeRegion;
   vtkSmartPointer<vtkDoubleArray> spanRegion;
   vtkSmartPointer<vtkCharArray>   typeRegion;

   int init(const idVertex numberOfVertices, Params params)
   {
      if (!params.segm)
         return 0;

      if (params.normalize) {
         normalizedIdVert = vtkSmartPointer<vtkIntArray>::New();
         normalizedIdVert->SetName("SegmentationId");
         normalizedIdVert->SetNumberOfComponents(1);
         normalizedIdVert->SetNumberOfTuples(numberOfVertices);
      } else {
         idVerts = vtkSmartPointer<vtkIntArray>::New();
         idVerts->SetName("SegmentationId");
         idVerts->SetNumberOfComponents(1);
         idVerts->SetNumberOfTuples(numberOfVertices);
      }

      sizeRegion = vtkSmartPointer<vtkIntArray>::New();
      sizeRegion->SetName("RegionSize");
      sizeRegion->SetNumberOfComponents(1);
      sizeRegion->SetNumberOfTuples(numberOfVertices);

      spanRegion = vtkSmartPointer<vtkDoubleArray>::New();
      spanRegion->SetName("RegionSpan");
      spanRegion->SetNumberOfComponents(1);
      spanRegion->SetNumberOfTuples(numberOfVertices);

      typeRegion = vtkSmartPointer<vtkCharArray>::New();
      typeRegion->SetName("RegionType");
      typeRegion->SetNumberOfComponents(1);
      typeRegion->SetNumberOfTuples(numberOfVertices);

// Check
#ifndef withKamikaze

      if (!idVerts && !normalizedIdVert) {
         cerr << "[ttkFTMTree] Error : vtkIntArray id allocation problem." << endl;
         return -2;
      }

      if (!sizeRegion) {
         cerr << "[ttkFTMTree] Error : vtkIntArray region size allocation problem." << endl;
         return -2;
      }

      if (!spanRegion) {
         cerr << "[ttkFTMTree] Error : vtkDoubleArray region span allocation problem." << endl;
         return -2;
      }

      if (!typeRegion) {
         cerr << "[ttkFTMTree] Error : vtkCharArray region type allocation problem." << endl;
         return -2;
      }

#endif
      return 0;
   }

   void fillArray(const idSuperArc arcId, Triangulation* triangulation, FTMTree_MT* tree,
                  Params params)
   {
      if (!params.segm)
         return;

      auto getNodeType = [&](const idNode nodeId) {
         Node* node = tree->getNode(nodeId);
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
      };

      SuperArc* arc = tree->getSuperArc(arcId);

      const int      upNodeId   = arc->getUpNodeId();
      const Node*    upNode     = tree->getNode(upNodeId);
      const int      upVertexId = upNode->getVertexId();
      const NodeType upNodeType = getNodeType(upNodeId);
      float          coordUp[3];
      triangulation->getVertexPoint(upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const int      downNodeId   = arc->getDownNodeId();
      const Node*    downNode     = tree->getNode(downNodeId);
      const int      downVertexId = downNode->getVertexId();
      const NodeType downNodeType = getNodeType(downNodeId);
      float          coordDown[3];
      triangulation->getVertexPoint(downVertexId, coordDown[0], coordDown[1], coordDown[2]);

      const int    regionSize = tree->getSuperArc(arcId)->getNumberOfRegularNodes();
      const double regionSpan = Geometry::distance(coordUp, coordDown);

      idSuperArc nid = arc->getNormalizedId();

      ArcType regionType;
      // RegionType
      if (upNodeType == NodeType::Local_minimum && downNodeType == NodeType::Local_maximum)
         regionType = ArcType::Min_arc;
      else if (upNodeType == NodeType::Local_minimum || downNodeType == NodeType::Local_minimum)
         regionType = ArcType::Min_arc;
      else if (upNodeType == NodeType::Local_maximum || downNodeType == NodeType::Local_maximum)
         regionType = ArcType::Max_arc;
      else if (upNodeType == NodeType::Saddle1 && downNodeType == NodeType::Saddle1)
         regionType = ArcType::Saddle1_arc;
      else if (upNodeType == NodeType::Saddle2 && downNodeType == NodeType::Saddle2)
         regionType = ArcType::Saddle2_arc;
      else
         regionType = ArcType::Saddle1_saddle2_arc;

      // fill extrema and regular verts of this arc

      // critical points
      if(params.normalize){
         normalizedIdVert->SetTuple1(upVertexId   , nid);
         normalizedIdVert->SetTuple1(downVertexId , nid);
      } else{
         idVerts->SetTuple1(upVertexId   , arcId);
         idVerts->SetTuple1(downVertexId , arcId);
      }

      sizeRegion->SetTuple1(upVertexId   , regionSize);
      sizeRegion->SetTuple1(downVertexId , regionSize);
      spanRegion->SetTuple1(upVertexId   , regionSpan);
      spanRegion->SetTuple1(downVertexId , regionSpan);
      typeRegion->SetTuple1(upVertexId   , static_cast<char>(regionType));
      typeRegion->SetTuple1(downVertexId , static_cast<char>(regionType));

      // regular nodes
      for (const idVertex vertexId : *arc) {
         if (params.normalize) {
            normalizedIdVert->SetTuple1(vertexId, nid);
         } else {
            idVerts->SetTuple1(vertexId, arcId);
         }
         sizeRegion->SetTuple1(vertexId, regionSize);
         spanRegion->SetTuple1(vertexId, regionSpan);
         typeRegion->SetTuple1(vertexId, static_cast<char>(regionType));
      }

   }

   void addArray(vtkPointData* pointData, Params params){
      if (!params.segm)
         return;

      if(params.normalize){
         pointData->AddArray(normalizedIdVert);
      } else {
         pointData->AddArray(idVerts);
      }

      pointData->AddArray(sizeRegion);
      pointData->AddArray(spanRegion);
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
       params_.treeType = (TreeType)type;
       Modified();
   }

   TreeType GetTreeType(void) const
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

   int setupTriangulation(vtkDataSet* input);
   int getScalars(vtkDataSet* input);
   int getOffsets(vtkDataSet* input);

   NodeType getNodeType(const Node* node);

   int getSkeletonNodes(FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(FTMTree_MT* tree,
                            idSuperArc arcId,
                            vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs,
                            ArcData& arcData);

   int addSampledSkeletonArc(FTMTree_MT* tree,
                             idSuperArc arcId,
                             const int samplingLevel,
                             vtkPoints* points,
                             vtkUnstructuredGrid* skeletonArcs,
                             ArcData& arcData);

   int addCompleteSkeletonArc(FTMTree_MT* tree,
                              idSuperArc arcId,
                              vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs,
                              ArcData& arcData);

   int getSkeletonArcs(FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonArcs);

   int getSegmentation(FTMTree_MT* tree, vtkDataSet* input, vtkDataSet* outputSegmentation);

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

   Params params_;

   Triangulation* triangulation_;
   FTMTree        ftmTree_;
   vtkDataArray*  inputScalars_;
   vtkIntArray*   offsets_;
   vtkDataArray*  inputOffsets_;
   bool           hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
