#pragma once

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include "FTRCommon.h"
#include "FTRDataTypes.h"
#include "Graph.h"

namespace ttk {
  namespace ftr {

    /// Vertex / Node / Arc data inherit from this
    /// master structure.
    struct ObjectData {
      inline void allocArray(vtkDataArray *const arr,
                             const char *fieldName,
                             size_t nbElmnt) {
        arr->SetName(fieldName);
        arr->SetNumberOfComponents(1);
        arr->SetNumberOfTuples(nbElmnt);

#ifndef TTK_ENABLE_KAMIKAZE
        if(!arr) {
          Debug dbg{};
          dbg.setDebugMsgPrefix("FTRGraph");
          dbg.printErr("unable to allocate " + std::string{fieldName}
                       + " the program will likely crash");
        }
#endif
      }
    };

    struct NodeData : public ObjectData {
      vtkNew<vtkIntArray> ids{};
      vtkNew<vtkIntArray> types{};
      vtkNew<vtkDoubleArray> scalars{};

      explicit NodeData(const ttk::ftr::idVertex nbNodes) {
        allocArray(ids, "VertexId", nbNodes);
        allocArray(types, "CriticalType", nbNodes);
        allocArray(scalars, "Scalar", nbNodes);
      }

      void addNode(const ttk::ftr::Graph &graph,
                   const ttk::ftr::idNode n,
                   const double scalar) {
        ids->SetTuple1(n, graph.getNode(n).getVertexIdentifier());
        types->SetTuple1(n, (double)graph.getNode(n).getType());
        scalars->SetTuple1(n, scalar);
      }

      void addArrays(vtkPointData *pointData,
                     ttk::ftr::Params ttkNotUsed(params)) {
        pointData->AddArray(ids);
        pointData->SetScalars(types);
        pointData->AddArray(scalars);
      }
    };

    struct ArcData : public ObjectData {
      vtkNew<vtkIntArray> ids{};
      vtkNew<vtkCharArray> reg{};
#ifndef NDEBUG
      vtkNew<vtkUnsignedCharArray> fromUp{};
#endif
      std::map<ttk::ftr::idVertex, vtkIdType> points;

      ArcData(const ttk::ftr::idSuperArc nbArcs) {
        allocArray(ids, "ArcId", nbArcs);
        allocArray(reg, ttk::MaskScalarFieldName, nbArcs * 2);
#ifndef NDEBUG
        allocArray(fromUp, "growUp", nbArcs);
#endif
      }

      void setPointInfo(const ttk::ftr::Graph &ttkNotUsed(graph),
                        const ttk::ftr::idSuperArc ttkNotUsed(a),
                        const vtkIdType skeletonVert,
                        bool r = false) {
        reg->SetTuple1(skeletonVert, r);
      }

      void setArcInfo(const ttk::ftr::Graph &graph,
                      const ttk::ftr::idSuperArc a,
                      const vtkIdType skeletonCell) {
        ids->SetTuple1(skeletonCell, a);
#ifndef NDEBUG
        fromUp->SetTuple1(skeletonCell, graph.getArc(a).getFromUp());
#else
        TTK_FORCE_USE(graph);
#endif
      }

      void addArrays(vtkUnstructuredGrid *arcs,
                     ttk::ftr::Params ttkNotUsed(params)) {
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
      vtkNew<vtkIntArray> ids{};
      vtkNew<vtkIntArray> regionType{};
#ifdef TTK_ENABLE_FTR_VERT_STATS
      vtkNew<vtkIntArray> touch{};
      vtkNew<vtkIntArray> arcActif{};
      vtkNew<vtkIntArray> taskActif{};
#endif

      explicit VertData(const ttk::ftr::idVertex nbVertices) {
        allocArray(ids, "ArcId", nbVertices);
        allocArray(regionType, "RegionType", nbVertices);
#ifdef TTK_ENABLE_FTR_VERT_STATS
        allocArray(touch, "Visit", nbVertices);
        allocArray(arcActif, "Arc active", nbVertices);
#endif
      }

      void setVertexInfo(const ttk::ftr::Graph &graph,
                         const ttk::ftr::idVertex v) {

        if(!graph.isVisited(v)) {
          // Problem, we should have visited all vertices
          // return to avoid crash
          return;
        }

        const ttk::ftr::idSuperArc curArcId = graph.getArcId(v);
        ids->SetTuple1(v, curArcId);

        int downNodeType
          = (int)graph.getNode(graph.getArc(curArcId).getDownNodeId())
              .getType();
        regionType->SetTuple1(v, downNodeType);

#ifdef TTK_ENABLE_FTR_VERT_STATS
        touch->SetTuple1(v, graph.getNbTouch(v));
        arcActif->SetTuple1(v, graph.getNbArcActive(v));
#endif
      }

      void addArrays(vtkDataSet *segmentation,
                     ttk::ftr::Params ttkNotUsed(params)) {
        segmentation->GetPointData()->AddArray(ids);
        segmentation->GetPointData()->SetActiveScalars(ids->GetName());
        segmentation->GetPointData()->AddArray(regionType);
#ifdef TTK_ENABLE_FTR_VERT_STATS
        segmentation->GetPointData()->AddArray(touch);
        segmentation->GetPointData()->AddArray(arcActif);
#endif
      }
    };
  }; // namespace ftr
}; // namespace ttk
