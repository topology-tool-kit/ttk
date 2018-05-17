#pragma once

#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkDataSet.h>

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

   explicit NodeData(const ttk::ftr::idVertex nbNodes)
   {
      ids = allocArray<vtkIntArray>("NodeId", nbNodes);
   }

   void addNode(const ttk::ftr::Graph& graph, const ttk::ftr::idNode n)
   {
      ids->SetTuple1(n, graph.getNode(n).getVertexIdentifier());
   }

   void addArrays(vtkPointData* pointData, ttk::ftr::Params params)
   {
      pointData->AddArray(ids);
   }
};

struct ArcData : public ObjectData {
   vtkSmartPointer<vtkIntArray> ids;
   std::map<ttk::ftr::idVertex, vtkIdType> points;

   explicit ArcData(const ttk::ftr::idSuperArc nbArcs)
   {
      ids = allocArray<vtkIntArray>("ArcId", nbArcs);
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
      segmentation->GetPointData()->AddArray(ids);
   }
};
