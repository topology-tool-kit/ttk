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
   vtkSmartPointer<vtkIntArray> sizeArcs;
#ifdef withStatsTime
   vtkSmartPointer<vtkFloatArray> startArcs;
   vtkSmartPointer<vtkFloatArray> endArcs;
   vtkSmartPointer<vtkFloatArray> timeArcs;
   vtkSmartPointer<vtkIntArray>   origArcs;
   vtkSmartPointer<vtkIntArray>   tasksArcs;
#endif

   int init()
   {
      sizeArcs = vtkSmartPointer<vtkIntArray>::New();
      sizeArcs->SetName("Size");

#ifdef withStatsTime
      startArcs = vtkSmartPointer<vtkFloatArray>::New();
      startArcs->SetName("Start");

      endArcs = vtkSmartPointer<vtkFloatArray>::New();
      endArcs->SetName("End");

      timeArcs = vtkSmartPointer<vtkFloatArray>::New();
      timeArcs->SetName("Time");

      origArcs = vtkSmartPointer<vtkIntArray>::New();
      origArcs->SetName("Origin");

      tasksArcs = vtkSmartPointer<vtkIntArray>::New();
      tasksArcs->SetName("Tasks");
#endif

// Check
#ifndef withKamikaze

      if (!sizeArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray size allocation problem." << endl;
         return -1;
      }
# ifdef withStatsTime
      if (!startArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray start allocation problem." << endl;
         return -2;
      }

      if (!endArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray end allocation problem." << endl;
         return -2;
      }

      if (!timeArcs) {
         cerr << "[ttkFTMTree] Error : vtkFloatArray time allocation problem." << endl;
         return -2;
      }

      if (!origArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray origin allocation problem." << endl;
         return -3;
      }

      if (!tasksArcs) {
         cerr << "[ttkFTMTree] Error : vtkIntArray tasks allocation problem." << endl;
         return -3;
      }
# endif
#endif
      return 0;
   }

   void addArray(vtkUnstructuredGrid *skeletonArcs){
      skeletonArcs->GetCellData()->AddArray(sizeArcs);
#ifdef withStatsTime
      skeletonArcs->GetCellData()->AddArray(startArcs);
      skeletonArcs->GetCellData()->AddArray(endArcs);
      skeletonArcs->GetCellData()->AddArray(timeArcs);
      skeletonArcs->GetCellData()->AddArray(origArcs);
      skeletonArcs->GetCellData()->AddArray(tasksArcs);
#endif
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

   vtkSetMacro(treeType_, int);
   vtkGetMacro(treeType_, int);

   vtkSetMacro(withSegmentation_, bool);
   vtkGetMacro(withSegmentation_, bool);

   vtkSetMacro(ScalarFieldId, int);
   vtkGetMacro(ScalarFieldId, int);

   vtkSetMacro(OffsetFieldId, int);
   vtkGetMacro(OffsetFieldId, int);

   vtkSetMacro(SuperArcSamplingLevel, int);
   vtkGetMacro(SuperArcSamplingLevel, int);

   int setupTriangulation(vtkDataSet* input);
   int getScalars(vtkDataSet* input);
   int getOffsets(vtkDataSet* input);

   NodeType getNodeType(const Node* node);

   int getSkeletonNodes(FTMTree_MT* tree, vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(FTMTree_MT* tree, SuperArc* arc, vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs);
   int addSampledSkeletonArc(FTMTree_MT* tree, SuperArc* arc, const int samplingLevel,
                             vtkPoints* points, vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);
   int addCompleteSkeletonArc(FTMTree_MT* tree, SuperArc* arc, vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs);
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

   int  treeType_;
   bool withSegmentation_;

   Triangulation* triangulation_;
   FTMTree        ftmTree_;
   vtkDataArray*  inputScalars_;
   vtkIntArray*   offsets_;
   vtkDataArray*  inputOffsets_;
   bool           hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
