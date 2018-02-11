#include "ttkFTRGraph.h"
#include "ttkFTRGraph.h"

// only used on the cpp
#include <vtkConnectivityFilter.h>
#include <vtkDataObject.h>
#include <vtkThreshold.h>

using namespace ttk::ftr;

vtkStandardNewMacro(ttkFTRGraph);

int ttkFTRGraph::FillInputPortInformation(int port, vtkInformation* info)
{
   if (port == 0)
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
   return 1;
}

int ttkFTRGraph::FillOutputPortInformation(int port, vtkInformation* info)
{
   switch (port) {
      case 0:
      case 1:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
         break;

      case 2:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
         break;
   }

   return 1;
}

int ttkFTRGraph::addCompleteSkeletonArc(const ttk::ftr::idSuperArc arcId, const int cc,
                                        vtkPoints* points, vtkUnstructuredGrid* skeletonArcs,
                                        ArcData& arcData)
{
   return 0;
}

int ttkFTRGraph::addDirectSkeletonArc(const idSuperArc arcId, const int cc, vtkPoints* points,
                                      vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   return 0;
}

int ttkFTRGraph::addSampledSkeletonArc(const idSuperArc arcId, const int cc, vtkPoints* points,
                                       vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   return 0;
}

int ttkFTRGraph::doIt(std::vector<vtkDataSet*>& inputs, std::vector<vtkDataSet*>& outputs)
{
   ttk::Memory m;

   // Input
   mesh_ = inputs[0];

   // Output
   vtkUnstructuredGrid* outputSkeletonNodes = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
   vtkUnstructuredGrid* outputSkeletonArcs  = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
   vtkDataSet*          outputSegmentation  = outputs[2];

   // Skeleton
   Graph graph;

   if (setupTriangulation()) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTRGraph] Error : wrong triangulation." << endl;
      return -1;
#endif
   }

   // Scalar field related
   getScalars();
   getOffsets();

   if (debugLevel_) {
      cout << "Launch on field : " << ScalarField << endl;
   }

   // compute graph
   switch (inputScalars_->GetDataType()) {
      vtkTemplateMacro({
         ttk::ftr::FTRGraph<VTK_TT> ftrGraph_;
         // common parameters
         ftrGraph_.setDebugLevel(debugLevel_);
         ftrGraph_.setThreadNumber(threadNumber_);
         ftrGraph_.setupTriangulation(triangulation_);
         // reeb graph parameters
         ftrGraph_.setScalars(inputScalars_->GetVoidPointer(0));
         ftrGraph_.setVertexSoSoffsets(offsets_.data());
         // build
         ftrGraph_.build();
         // get output
         graph = std::move(ftrGraph_.extractOutputGraph());
      });
   }

   UpdateProgress(0.50);

   // Construct output
   if (getSkeletonNodes(graph, outputSkeletonNodes)) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTRGraph] Error : wrong properties on skeleton nodes." << endl;
      return -7;
#endif
   }

   if (getSkeletonArcs(graph, outputSkeletonArcs)) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTRGraph] Error : wrong properties on skeleton arcs." << endl;
      return -8;
#endif
   }

   if (GetWithSegmentation()) {
      outputSegmentation->ShallowCopy(inputs[0]);
      if (getSegmentation(graph, outputSegmentation)) {
#ifndef TTK_ENABLE_KAMIKAZE
         cerr << "[ttkFTRGraph] Error : wrong properties on segmentation." << endl;
         return -9;
#endif
      }
   }

   UpdateProgress(1);

   {
      std::stringstream msg;
      msg << "[ttkFTRGraph] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}

int ttkFTRGraph::getOffsets()
{
   vtkDataArray* inputOffsets;
   inputOffsets = mesh_->GetPointData()->GetArray(OffsetFieldId);
   if (inputOffsets) {
      InputOffsetScalarFieldName = inputOffsets->GetName();
      UseInputOffsetScalarField  = true;
   }

   const idVertex numberOfVertices = mesh_->GetNumberOfPoints();

   if (UseInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
      inputOffsets = mesh_->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
      offsets_.resize(numberOfVertices);
      for (int i = 0; i < numberOfVertices; i++) {
         offsets_[i] = inputOffsets->GetTuple1(i);
      }
   } else {
      if (hasUpdatedMesh_ and offsets_.size()) {
         // don't keep an out-dated offset array
         offsets_.clear();
      }

      if (offsets_.empty()) {
         offsets_.resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static, numberOfVertices / threadNumber_)
#endif
         for (int i = 0; i < numberOfVertices; i++) {
            offsets_[i] = i;
         }
      }
   }

#ifndef TTK_ENABLE_KAMIKAZE
   if (offsets_.empty()) {
      cerr << "[ttkFTRGraph] Error : wrong input offset scalar field " << endl;
      return -1;
   }
#endif

   return 0;
}

int ttkFTRGraph::getScalars()
{
   vtkPointData* pointData = mesh_->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
   if (!pointData) {
      cerr << "[ttkFTRGraph] Error : input has no point data." << endl;
      return -1;
   }
#endif

   if (ScalarField.length()) {
      inputScalars_ = pointData->GetArray(ScalarField.data());
   } else {
      inputScalars_ = pointData->GetArray(ScalarFieldId);
      if (inputScalars_)
         ScalarField = inputScalars_->GetName();
   }

#ifndef TTK_ENABLE_KAMIKAZE
   if (!inputScalars_) {
      cerr << "[ttkFTRGraph] Error : input scalar field pointer is null." << endl;
      return -3;
   }
#endif

   return 0;
}

int ttkFTRGraph::getSegmentation(const ttk::ftr::Graph& graph, vtkDataSet* outputSegmentation)
{
   const idVertex numberOfVertices = mesh_->GetNumberOfPoints();
   VertData vertData(numberOfVertices);

   // TODO parallel
   for (idVertex v = 0; v < numberOfVertices; ++v) {
      vertData.setVertexInfo(graph, v);
   }

   vertData.addArrays(outputSegmentation, params_);

   return 0;
}

int ttkFTRGraph::getSkeletonArcs(const ttk::ftr::Graph& graph, vtkUnstructuredGrid* outputSkeletonArcs)
{
   return 0;
}

int ttkFTRGraph::getSkeletonNodes(const Graph& graph, vtkUnstructuredGrid* outputSkeletonNodes)
{
   const idNode nbNodes = graph.getNumberOfNodes();

   NodeData nodeData(nbNodes);
   vtkSmartPointer<vtkUnstructuredGrid> nodes  = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints>           points = vtkSmartPointer<vtkPoints>::New();

   for (idNode i = 0; i < nbNodes; i++) {
      const ttk::ftr::idVertex vertId = graph.getNode(i).getVertexIdentifier();
      float point[3];
      triangulation_->getVertexPoint(vertId, point[0], point[1], point[2]);
      points->InsertNextPoint(point);

      nodeData.addNode(graph, i);
   }

   nodes->SetPoints(points);
   nodeData.addArrays(nodes->GetPointData(), params_);

   outputSkeletonNodes->ShallowCopy(nodes);
   return 0;
}

int ttkFTRGraph::setupTriangulation()
{
   triangulation_ = ttkTriangulation::getTriangulation(mesh_);
#ifndef TTK_ENABLE_KAMIKAZE
   if (!triangulation_) {
      cerr << "[ttkFTRGraph] Error : ttkTriangulation::getTriangulation() is null." << endl;
      return -1;
   }
#endif

   triangulation_->setWrapper(this);
   hasUpdatedMesh_ = ttkTriangulation::hasChangedConnectivity(triangulation_, mesh_, this);

#ifndef TTK_ENABLE_KAMIKAZE
   if (triangulation_->isEmpty()) {
      cerr << "[ttkFTRGraph] Error : ttkTriangulation on connected component allocation problem."
           << endl;
      return -1;
   }
#endif
   return 0;
}

// protected

ttkFTRGraph::ttkFTRGraph()
    : ScalarField{},
      UseInputOffsetScalarField{},
      InputOffsetScalarFieldName{},
      ScalarFieldId{},
      OffsetFieldId{-1},
      params_{},
      mesh_{},
      triangulation_{},
      inputScalars_{},
      offsets_{},
      hasUpdatedMesh_{}
{
   SetNumberOfInputPorts(1);
   SetNumberOfOutputPorts(3);
}

ttkFTRGraph::~ttkFTRGraph()
{
}

void ttkFTRGraph::identify(vtkDataSet* ds) const
{
   vtkSmartPointer<vtkIntArray> identifiers = vtkSmartPointer<vtkIntArray>::New();
   const vtkIdType              nbPoints    = ds->GetNumberOfPoints();
   identifiers->SetName("VertexIdentifier");
   identifiers->SetNumberOfComponents(1);
   identifiers->SetNumberOfTuples(nbPoints);

   for (int i = 0; i < nbPoints; i++) {
      identifiers->SetTuple1(i, i);
   }

   ds->GetPointData()->AddArray(identifiers);
}
