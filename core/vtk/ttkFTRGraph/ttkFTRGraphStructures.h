#pragma once

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "DataTypesFTR.h"
#include "FTRCommon.h"
#include "Graph.h"

/// Vertex / Node / Arc data inherit from this
/// master structure.
struct ObjectData {
   template <typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> allocArray(const char* fieldName, size_t nbElmnt)
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
};

struct NodeData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkIntArray> types;
   vtkSmartPointer<vtkDoubleArray> scalars;

   explicit NodeData(const ttk::ftr::idVertex nbNodes)
   {
      ids     = allocArray<vtkIntArray>("VertexId", nbNodes);
      types   = allocArray<vtkIntArray>("NodeType", nbNodes);
      scalars = allocArray<vtkDoubleArray>("Scalar", nbNodes);
   }

   void addNode(const ttk::ftr::Graph& graph, const ttk::ftr::idNode n, const double scalar)
   {
      ids->SetTuple1(n, graph.getNode(n).getVertexIdentifier());
      types->SetTuple1(n, (double)graph.getNode(n).getType());
      scalars->SetTuple1(n, scalar);
   }

   void addArrays(vtkPointData* pointData, ttk::ftr::Params params)
   {
      pointData->AddArray(ids);
      pointData->SetScalars(types);
      pointData->AddArray(scalars);
   }
};

struct ArcData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkCharArray> reg;
#ifndef NDEBUG
   vtkSmartPointer<vtkUnsignedCharArray> fromUp;
#endif
   std::map<ttk::ftr::idVertex, vtkIdType> points; 
                                                   
   explicit ArcData(const ttk::ftr::idSuperArc nbArcs)
   {                                               
      ids = allocArray<vtkIntArray>("ArcId", nbArcs);
      reg = allocArray<vtkCharArray>("RegularMask", nbArcs+1);
#ifndef NDEBUG                                     
      fromUp = allocArray<vtkUnsignedCharArray>("growUp", nbArcs);
#endif
   }

   void setPointInfo(const ttk::ftr::Graph& graph, const ttk::ftr::idSuperArc a,
                     const vtkIdType skeletonVert, bool r = false)
   {
      reg->SetTuple1(skeletonVert, r);
   }

   void setArcInfo(const ttk::ftr::Graph& graph, const ttk::ftr::idSuperArc a,
                   const vtkIdType skeletonCell)
   {
      ids->SetTuple1(skeletonCell, a);
#ifndef NDEBUG
      fromUp->SetTuple1(skeletonCell, graph.getArc(a).getFromUp());
#endif
   }

   void addArrays(vtkUnstructuredGrid* arcs, ttk::ftr::Params params)
   {
      // original size may be too large
      ids->SetNumberOfTuples(arcs->GetNumberOfCells());
      arcs->GetCellData()->SetScalars(ids);
      reg->SetNumberOfTuples(arcs->GetNumberOfPoints());
      arcs->GetPointData()->AddArray(reg);
#ifndef NDEBUG
      fromUp->SetNumberOfTuples(arcs->GetNumberOfCells());
      arcs->GetCellData()->AddArray(fromUp);
#endif
   }
};

struct VertData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;
#ifdef TTK_ENABLE_FTR_VERT_STATS
   vtkSmartPointer<vtkIntArray> touch;
   vtkSmartPointer<vtkIntArray> arcActif;
   vtkSmartPointer<vtkIntArray> taskActif;
#endif

   explicit VertData(const ttk::ftr::idVertex nbVertices)
   {
      ids = allocArray<vtkIntArray>("ArcId", nbVertices);
#ifdef TTK_ENABLE_FTR_VERT_STATS
      touch     = allocArray<vtkIntArray>("Visit", nbVertices);
      arcActif  = allocArray<vtkIntArray>("Arc active", nbVertices);
#endif
   }

   void setVertexInfo(const ttk::ftr::Graph& graph, const ttk::ftr::idVertex v)
   {
      if (graph.isVisited(v)) {
         ids->SetTuple1(v, graph.getArcId(v));
      } else {
         ids->SetTuple1(v, -1);
      }

#ifdef TTK_ENABLE_FTR_VERT_STATS
      touch->SetTuple1(v, graph.getNbTouch(v));
      arcActif->SetTuple1(v, graph.getNbArcActive(v));
#endif
   }

   void addArrays(vtkDataSet* segmentation, ttk::ftr::Params params)
   {
      segmentation->GetPointData()->SetScalars(ids);
#ifdef TTK_ENABLE_FTR_VERT_STATS
      segmentation->GetPointData()->AddArray(touch);
      segmentation->GetPointData()->AddArray(arcActif);
#endif
   }
};
