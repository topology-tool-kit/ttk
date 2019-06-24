#include <ttkManifoldCheck.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkManifoldCheck)

  int ttkManifoldCheck::doIt(vector<vtkDataSet *> &inputs,
                             vector<vtkDataSet *> &outputs) {

  Memory mem;

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  triangulation->setWrapper(this);
  triangulation->preprocessVertexTriangles();
  manifoldCheck_.setupTriangulation(triangulation);
  manifoldCheck_.setWrapper(this);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  manifoldCheck_.setVertexLinkComponentNumberVector(
    &vertexLinkComponentNumber_);
  manifoldCheck_.setEdgeLinkComponentNumberVector(&edgeLinkComponentNumber_);
  manifoldCheck_.setTriangleLinkComponentNumberVector(
    &triangleLinkComponentNumber_);
  manifoldCheck_.execute();

  {
    stringstream msg;
    msg << "[ttkManifoldCheck] Preparing VTK output..." << endl;
    dMsg(cout, msg.str(), Debug::timeMsg);
  }

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexPointArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexPointArray->SetName("VertexLinkComponentNumber");
  vertexPointArray->SetNumberOfTuples(output->GetNumberOfPoints());
  for(SimplexId i = 0; i < (SimplexId)vertexLinkComponentNumber_.size(); i++)
    vertexPointArray->SetTuple1(i, vertexLinkComponentNumber_[i]);
  output->GetPointData()->AddArray(vertexPointArray);

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexCellArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexCellArray->SetName("VertexLinkComponentNumber");
  vertexCellArray->SetNumberOfTuples(output->GetNumberOfCells());

  for(SimplexId i = 0; i < output->GetNumberOfCells(); i++) {
    vtkCell *c = output->GetCell(i);
    SimplexId cellMax = -1;
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      SimplexId vertexId = c->GetPointId(j);
      if((!j) || (vertexLinkComponentNumber_[vertexId] > cellMax)) {
        cellMax = vertexLinkComponentNumber_[vertexId];
      }
    }

    vertexCellArray->SetTuple1(i, cellMax);
  }
  output->GetCellData()->AddArray(vertexCellArray);

  // edges
  vtkSmartPointer<ttkSimplexIdTypeArray> edgePointArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  edgePointArray->SetName("EdgeLinkComponentNumber");
  edgePointArray->SetNumberOfTuples(output->GetNumberOfPoints());
  for(SimplexId i = 0; i < edgePointArray->GetNumberOfTuples(); i++) {
    edgePointArray->SetTuple1(i, 0);
  }

  vtkSmartPointer<ttkSimplexIdTypeArray> edgeCellArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  edgeCellArray->SetName("EdgeLinkComponentNumber");
  edgeCellArray->SetNumberOfTuples(output->GetNumberOfCells());
  for(SimplexId i = 0; i < edgeCellArray->GetNumberOfTuples(); i++) {
    edgeCellArray->SetTuple1(i, 0);
  }

  if(edgeLinkComponentNumber_.size()) {

    for(SimplexId i = 0; i < (SimplexId)edgeLinkComponentNumber_.size(); i++) {

      SimplexId vertexId0 = -1, vertexId1 = -1;
      triangulation->getEdgeVertex(i, 0, vertexId0);
      triangulation->getEdgeVertex(i, 1, vertexId1);

      SimplexId vertexMax0 = edgePointArray->GetTuple1(vertexId0);
      SimplexId vertexMax1 = edgePointArray->GetTuple1(vertexId1);

      if(edgeLinkComponentNumber_[i] > vertexMax0)
        edgePointArray->SetTuple1(vertexId0, edgeLinkComponentNumber_[i]);
      if(edgeLinkComponentNumber_[i] > vertexMax1)
        edgePointArray->SetTuple1(vertexId1, edgeLinkComponentNumber_[i]);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < output->GetNumberOfCells(); i++) {
      vtkSmartPointer<vtkGenericCell> c
        = vtkSmartPointer<vtkGenericCell>::New();
      output->GetCell(i, c);
      SimplexId cellMax = -1;
      for(int j = 0; j < c->GetNumberOfPoints(); j++) {
        SimplexId vertexId0 = c->GetPointId(j);
        SimplexId vertexId1 = -1;
        for(int k = 0; k < c->GetNumberOfPoints(); k++) {
          if(k != j) {
            vertexId1 = c->GetPointId(k);

            // check if (vertexId0 - vertexId1) is indeed an edge in the
            // triangulation
            SimplexId edgeNumber
              = triangulation->getVertexEdgeNumber(vertexId0);
            for(SimplexId l = 0; l < edgeNumber; l++) {
              SimplexId edgeId = -1;
              triangulation->getVertexEdge(vertexId0, l, edgeId);

              SimplexId vertexIdA = -1, vertexIdB = -1;
              triangulation->getEdgeVertex(edgeId, 0, vertexIdA);
              triangulation->getEdgeVertex(edgeId, 1, vertexIdB);

              if(((vertexId0 == vertexIdA) && (vertexId1 == vertexIdB))
                 || ((vertexId1 == vertexIdA) && (vertexId0 == vertexIdB))) {

                // (vertexId0 - vertexId1) is indeed an edge in the
                // triangulation
                if(edgeLinkComponentNumber_[edgeId] > cellMax)
                  cellMax = edgeLinkComponentNumber_[edgeId];
              }
            }
          }
        }
      }
      edgeCellArray->SetTuple1(i, cellMax);
    }
  }
  output->GetPointData()->AddArray(edgePointArray);
  output->GetCellData()->AddArray(edgeCellArray);

  // triangles
  vtkSmartPointer<ttkSimplexIdTypeArray> trianglePointArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  trianglePointArray->SetName("TriangleLinkComponentNumber");
  trianglePointArray->SetNumberOfTuples(output->GetNumberOfPoints());
  for(SimplexId i = 0; i < trianglePointArray->GetNumberOfTuples(); i++) {
    trianglePointArray->SetTuple1(i, 0);
  }

  vtkSmartPointer<ttkSimplexIdTypeArray> triangleCellArray
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  triangleCellArray->SetName("TriangleLinkComponentNumber");
  triangleCellArray->SetNumberOfTuples(output->GetNumberOfCells());
  for(SimplexId i = 0; i < triangleCellArray->GetNumberOfTuples(); i++) {
    triangleCellArray->SetTuple1(i, 0);
  }

  if(triangleLinkComponentNumber_.size()) {

    for(SimplexId i = 0; i < (SimplexId)triangleLinkComponentNumber_.size();
        i++) {

      SimplexId vertexId0 = -1, vertexId1 = -1, vertexId2 = -1;
      triangulation->getTriangleVertex(i, 0, vertexId0);
      triangulation->getTriangleVertex(i, 1, vertexId1);
      triangulation->getTriangleVertex(i, 2, vertexId2);

      SimplexId vertexMax0 = trianglePointArray->GetTuple1(vertexId0);
      SimplexId vertexMax1 = trianglePointArray->GetTuple1(vertexId1);
      SimplexId vertexMax2 = trianglePointArray->GetTuple1(vertexId2);

      if(triangleLinkComponentNumber_[i] > vertexMax0)
        trianglePointArray->SetTuple1(
          vertexId0, triangleLinkComponentNumber_[i]);
      if(triangleLinkComponentNumber_[i] > vertexMax1)
        trianglePointArray->SetTuple1(
          vertexId1, triangleLinkComponentNumber_[i]);
      if(triangleLinkComponentNumber_[i] > vertexMax2)
        trianglePointArray->SetTuple1(
          vertexId2, triangleLinkComponentNumber_[i]);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < output->GetNumberOfCells(); i++) {
      vtkSmartPointer<vtkGenericCell> c
        = vtkSmartPointer<vtkGenericCell>::New();
      output->GetCell(i, c);

      SimplexId cellMax = -1;
      for(int j = 0; j < c->GetNumberOfPoints(); j++) {
        SimplexId vertexId0 = c->GetPointId(j);
        SimplexId vertexId1 = -1;
        SimplexId vertexId2 = -1;

        for(int k = 0; k < c->GetNumberOfPoints(); k++) {
          if(k != j) {
            vertexId1 = c->GetPointId(k);

            for(int l = 0; l < c->GetNumberOfPoints(); l++) {
              if((l != j) && (l != k)) {
                vertexId2 = c->GetPointId(l);

                // check if (vertexId0, vertexId1, vertexId2) is indeed a
                // triangle in the triangulation
                SimplexId triangleNumber
                  = triangulation->getVertexTriangleNumber(vertexId0);
                for(SimplexId m = 0; m < triangleNumber; m++) {
                  SimplexId triangleId = -1;
                  triangulation->getVertexTriangle(vertexId0, m, triangleId);

                  SimplexId vertexIdA = -1, vertexIdB = -1, vertexIdC = -1;
                  triangulation->getTriangleVertex(triangleId, 0, vertexIdA);
                  triangulation->getTriangleVertex(triangleId, 1, vertexIdB);
                  triangulation->getTriangleVertex(triangleId, 2, vertexIdC);

                  if(((vertexId0 == vertexIdA) && (vertexId1 == vertexIdB)
                      && (vertexId2 == vertexIdC))
                     /// ABC
                     || ((vertexId0 == vertexIdA) && (vertexId1 == vertexIdC)
                         && (vertexId2 == vertexIdB))
                     // ACB
                     || ((vertexId0 == vertexIdB) && (vertexId1 == vertexIdA)
                         && (vertexId2 == vertexIdC))
                     /// BAC
                     || ((vertexId0 == vertexIdB) && (vertexId1 == vertexIdC)
                         && (vertexId2 == vertexIdA))
                     // BCA
                     || ((vertexId0 == vertexIdC) && (vertexId1 == vertexIdA)
                         && (vertexId2 == vertexIdB))
                     /// CAB
                     || ((vertexId0 == vertexIdC) && (vertexId1 == vertexIdB)
                         && (vertexId2 == vertexIdA))) {
                    // CBA

                    // (vertexId0, vertexId1, vertexId2) is indeed a
                    // triangle in the triangulation
                    if(triangleLinkComponentNumber_[triangleId] > cellMax) {
                      cellMax = triangleLinkComponentNumber_[triangleId];
                    }
                  }
                }
              }
            }
          }
        }
      }
      triangleCellArray->SetTuple1(i, cellMax);
    }
  }
  output->GetPointData()->AddArray(trianglePointArray);
  output->GetCellData()->AddArray(triangleCellArray);

  {
    stringstream msg;
    msg << "[ttkManifoldCheck] Memory usage: " << mem.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
