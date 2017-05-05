/// \sa ttk::TaskedTree
#ifndef _VTK_CONTOURFORESTS_H
#define _VTK_CONTOURFORESTS_H

// ttk code includes
#include <TaskedTree.h>
#include <ttkWrapper.h>

// VTK includes
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

class VTKFILTERSCORE_EXPORT vtkTaskedTree : public vtkDataSetAlgorithm, public Wrapper
{
  public:
   static vtkTaskedTree* New();

   vtkTypeMacro(vtkTaskedTree, vtkDataSetAlgorithm);

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

   vtkSetMacro(lessPartition_, int);
   vtkGetMacro(lessPartition_, int);

   vtkSetMacro(partitionNumber_, int);
   vtkGetMacro(partitionNumber_, int);

   vtkSetMacro(simplificationMethod_, int);
   vtkGetMacro(simplificationMethod_, int);

   vtkSetMacro(simplificationThreshold_, double);
   vtkGetMacro(simplificationThreshold_, double);

   vtkSetMacro(useThresholdNormalization_, int);
   vtkGetMacro(useThresholdNormalization_, int);

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

   int getSkeletonNodes(MergeTree* tree, vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(MergeTree* tree, SuperArc* arc, vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs);
   int addSampledSkeletonArc(MergeTree* tree, SuperArc* arc, const int samplingLevel,
                             vtkPoints* points, vtkUnstructuredGrid* skeletonArcs);
   int addCompleteSkeletonArc(MergeTree* tree, SuperArc* arc, vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs);
   int getSkeletonArcs(MergeTree* tree, vtkUnstructuredGrid* outputSkeletonArcs);

   int getSegmentation(MergeTree* tree, vtkDataSet* input, vtkDataSet* outputSegmentation);

  protected:
   vtkTaskedTree();
   ~vtkTaskedTree();

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

   int    treeType_;
   int    lessPartition_;
   int    partitionNumber_;
   int    simplificationMethod_;
   double simplificationThreshold_;
   bool   useThresholdNormalization_;

   Triangulation* triangulation_;
   TaskedTree     contourForests_;
   vtkDataArray*  inputScalars_;
   vtkIntArray*   offsets_;
   vtkDataArray*  inputOffsets_;
   bool           hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
