/// \sa ttk::ftr::FTMTree
#ifndef _VTK_CONTOURFORESTS_H
#define _VTK_CONTOURFORESTS_H

// ttk code includes
#include "FTRGraph.h"
#include "Graph.h"
#include "ttkFTRGraphStructures.h"
#include "ttkWrapper.h"

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
// Data array
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

// Unused ? (compile without these on my computer)
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkObjectFactory.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFTRGraph : public vtkDataSetAlgorithm, public ttk::Wrapper
#else
class ttkFTRGraph : public vtkDataSetAlgorithm, public ttk::Wrapper
#endif
{
  private:
   std::string ScalarField;
   bool        UseInputOffsetScalarField;
   std::string InputOffsetScalarFieldName;
   int         ScalarFieldId;
   int         OffsetFieldId;

   ttk::ftr::Params params_;

   vtkDataSet*                     mesh_;
   ttk::Triangulation*             triangulation_;
   vtkDataArray*                   inputScalars_;
   std::vector<ttk::ftr::idVertex> offsets_;

   bool hasUpdatedMesh_;
   bool Dummy;

  public:
   static ttkFTRGraph* New();

   vtkTypeMacro(ttkFTRGraph, vtkDataSetAlgorithm);

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

   vtkSetMacro(ScalarField, std::string);
   vtkGetMacro(ScalarField, std::string);

   vtkSetMacro(UseInputOffsetScalarField, int);
   vtkGetMacro(UseInputOffsetScalarField, int);

   vtkSetMacro(InputOffsetScalarFieldName, std::string);
   vtkGetMacro(InputOffsetScalarFieldName, std::string);

   vtkSetMacro(ScalarFieldId, int);
   vtkGetMacro(ScalarFieldId, int);

   vtkSetMacro(OffsetFieldId, int);
   vtkGetMacro(OffsetFieldId, int);

   vtkSetMacro(Dummy, bool);

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

   int setupTriangulation();
   int getScalars();
   int getOffsets();

   int getSkeletonNodes(const ttk::ftr::Graph& graph, vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(const ttk::ftr::idSuperArc arcId, const int cc, vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int addSampledSkeletonArc(const ttk::ftr::idSuperArc arcId, const int cc, vtkPoints* points,
                             vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int addCompleteSkeletonArc(const ttk::ftr::idSuperArc arcId, const int cc, vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int getSkeletonArcs(const ttk::ftr::Graph& graph, vtkUnstructuredGrid* outputSkeletonArcs);

   int getSegmentation(const ttk::ftr::Graph& graph, vtkDataSet* outputSegmentation);

  protected:
   ttkFTRGraph();
   ~ttkFTRGraph();

   TTK_SETUP();

   void identify(vtkDataSet* ds) const;

   virtual int FillInputPortInformation(int port, vtkInformation* info) override;
   virtual int FillOutputPortInformation(int port, vtkInformation* info) override;
};

#endif  // _VTK_CONTOURFORESTS_H
