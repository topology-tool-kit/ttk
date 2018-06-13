#pragma once

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
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

   explicit NodeData(const ttk::ftr::idVertex nbNodes)
   {
      ids   = allocArray<vtkIntArray>("VertexId", nbNodes);
      types = allocArray<vtkIntArray>("NodeType", nbNodes);
   }

   void addNode(const ttk::ftr::Graph& graph, const ttk::ftr::idNode n)
   {
      ids->SetTuple1(n, graph.getNode(n).getVertexIdentifier());
      types->SetTuple1(n, (double)graph.getNode(n).getType());
   }

   void addArrays(vtkPointData* pointData, ttk::ftr::Params params)
   {
      pointData->AddArray(ids);
      pointData->SetScalars(types);
   }
};

struct ArcData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;
#ifndef NDEBUG
   vtkSmartPointer<vtkUnsignedCharArray> fromUp;
#endif
   std::map<ttk::ftr::idVertex, vtkIdType> points;

   explicit ArcData(const ttk::ftr::idSuperArc nbArcs)
   {
      ids = allocArray<vtkIntArray>("ArcId", nbArcs);
#ifndef NDEBUG
      fromUp = allocArray<vtkUnsignedCharArray>("growUp", nbArcs);
#endif
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
       arcs->GetCellData()->SetScalars(ids);
#ifndef NDEBUG
      arcs->GetCellData()->AddArray(fromUp);
#endif
   }
};

struct VertData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;

   explicit VertData(const ttk::ftr::idVertex nbVertices)
   {
      ids = allocArray<vtkIntArray>("Segmentation", nbVertices);
   }

   void setVertexInfo(const ttk::ftr::Graph& graph, const ttk::ftr::idVertex v)
   {
      if (graph.isVisited(v)) {
         ids->SetTuple1(v, graph.getFirstArcId(v));
      } else {
         ids->SetTuple1(v, -1);
      }
   }

   void addArrays(vtkDataSet* segmentation, ttk::ftr::Params params)
   {
      segmentation->GetPointData()->SetScalars(ids);
   }
};
