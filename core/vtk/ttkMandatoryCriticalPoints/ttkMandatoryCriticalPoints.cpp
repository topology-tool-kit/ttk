#include <ttkMandatoryCriticalPoints.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMandatoryCriticalPoints)

  ttkMandatoryCriticalPoints::ttkMandatoryCriticalPoints() {
  SetLowerBoundFieldName("lowerBoundField");
  SetUpperBoundFieldName("upperBoundField");
  SetSimplificationThreshold(0.0);

  // init
  outputMandatoryMinimum_ = NULL;
  outputMandatoryJoinSaddle_ = NULL;
  outputMandatorySplitSaddle_ = NULL;
  outputMandatoryMaximum_ = NULL;

  mandatoryJoinTreePoints_ = NULL;
  mandatorySplitTreePoints_ = NULL;

  UseAllCores = true;

  computeMinimumOutput_ = false;
  computeJoinSaddleOutput_ = false;
  computeSplitSaddleOutput_ = false;
  computeMaximumOutput_ = false;

  mandatoryJoinTreePoints_ = NULL;
  mdtJoinTreePointType_ = NULL;
  mdtJoinTreePointLowInterval_ = NULL;
  mdtJoinTreePointUpInterval_ = NULL;
  mdtJoinTreePointComponentId_ = NULL;
  mandatorySplitTreePoints_ = NULL;
  mdtSplitTreePointType_ = NULL;
  mdtSplitTreePointLowInterval_ = NULL;
  mdtSplitTreePointUpInterval_ = NULL;
  mdtSplitTreePointComponentId_ = NULL;
  mdtJoinTreeEdgeSwitchable_ = NULL;
  mdtSplitTreeEdgeSwitchable_ = NULL;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(6);

  inputMTime_ = 0;
  computeAll_ = true;
  simplificationThreshold_ = 0.0;
  simplify_ = true;

  memoryUsage_ = 0.0;

  lowerBoundId = 0;
  upperBoundId = 1;

  outputAllMinimumComponents_ = true;
  outputMinimumComponentId_ = 0;
  outputAllJoinSaddleComponents_ = true;
  outputJoinSaddleComponentId_ = 0;
  outputAllSplitSaddleComponents_ = true;
  outputSplitSaddleComponentId_ = 0;
  outputAllMaximumComponents_ = true;
  outputMaximumComponentId_ = 0;

  computeMinimumOutput_ = true;
  computeJoinSaddleOutput_ = true;
  computeSplitSaddleOutput_ = true;
  computeMaximumOutput_ = true;

  triangulation_ = NULL;
}

ttkMandatoryCriticalPoints::~ttkMandatoryCriticalPoints() {
  if(outputMandatoryMinimum_)
    outputMandatoryMinimum_->Delete();
  if(outputMandatoryJoinSaddle_)
    outputMandatoryJoinSaddle_->Delete();
  if(outputMandatorySplitSaddle_)
    outputMandatorySplitSaddle_->Delete();
  if(outputMandatoryMaximum_)
    outputMandatoryMaximum_->Delete();
  if(mandatoryJoinTreePoints_)
    mandatoryJoinTreePoints_->Delete();
  if(mdtJoinTreePointType_)
    mdtJoinTreePointType_->Delete();
  if(mdtJoinTreePointLowInterval_)
    mdtJoinTreePointLowInterval_->Delete();
  if(mdtJoinTreePointUpInterval_)
    mdtJoinTreePointUpInterval_->Delete();
  if(mdtJoinTreePointComponentId_)
    mdtJoinTreePointComponentId_->Delete();
  if(mandatorySplitTreePoints_)
    mandatorySplitTreePoints_->Delete();
  if(mdtSplitTreePointType_)
    mdtSplitTreePointType_->Delete();
  if(mdtSplitTreePointLowInterval_)
    mdtSplitTreePointLowInterval_->Delete();
  if(mdtSplitTreePointUpInterval_)
    mdtSplitTreePointUpInterval_->Delete();
  if(mdtSplitTreePointComponentId_)
    mdtSplitTreePointComponentId_->Delete();
  if(mandatoryJoinTreeEdge_.size()) {
    for(size_t i = 0; i < mandatoryJoinTreeEdge_.size(); i++) {
      if(mandatoryJoinTreeEdge_[i])
        mandatoryJoinTreeEdge_[i]->Delete();
    }
  }
  if(mandatorySplitTreeEdge_.size()) {
    for(size_t i = 0; i < mandatorySplitTreeEdge_.size(); i++) {
      if(mandatorySplitTreeEdge_[i])
        mandatorySplitTreeEdge_[i]->Delete();
    }
  }
  if(mdtJoinTreeEdgeSwitchable_)
    mdtJoinTreeEdgeSwitchable_->Delete();
  if(mdtSplitTreeEdgeSwitchable_)
    mdtSplitTreeEdgeSwitchable_->Delete();
}

// struct XCoordCmp {
//   bool operator()(const pair<int,double> &left, const pair<int,double>
//   &right){
//     return (left.second < right.second);
//   }
// } xCoordCmp;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// void ttkMandatoryCriticalPoints::buildMandatorySplitTree(vtkUnstructuredGrid
// *outputSplitTree) {
//
//   // Clear the objects
//   if(mandatorySplitTreePoints_)
//     mandatorySplitTreePoints_->Delete();
//   if(mandatorySplitTreeEdge_.size()) {
//     for(int i=0 ; i<(int)mandatorySplitTreeEdge_.size() ; i++)
//       if(mandatorySplitTreeEdge_[i])
//         mandatorySplitTreeEdge_[i]->Delete();
//     mandatorySplitTreeEdge_.clear();
//   }
//   // Tree Informations
//   const Digraph &mandatorySplitTree =
// mandatoryCriticalPoints_.getMandatorySplitTree();
//   /* New Objects */
//   // Points
//   mandatorySplitTreePoints_ = vtkPoints::New();
//   // Edges
//   mandatorySplitTreeEdge_.resize(mandatorySplitTree.getNumberOfEdges(),NULL);
//   for(int i=0 ; i<(int)mandatorySplitTreeEdge_.size() ; i++)
//     mandatorySplitTreeEdge_[i] = vtkIdList::New();
//
//   cout << "[ttkMandatoryCriticalPoints] Number of points : ";
//   cout << mandatorySplitTree.getNumberOfPoints() << endl;
//   cout << "[ttkMandatoryCriticalPoints] Number of edges : ";
//   cout << mandatorySplitTreeEdge_.size() << endl;
//
//   // PointId to GraphPointId
//   vector<int> graphPointToVtkPoint(mandatorySplitTree.getNumberOfPoints(),
// -1);
//   // Level
//   vector<int> pointLevel(mandatorySplitTree.getNumberOfPoints(), -1);
//   // Find the root
//   int rootGraphPointId = -1;
//   for(int i=0 ; i<(mandatorySplitTree.getNumberOfPoints()) ; i++) {
//     if(!mandatorySplitTree.getPoint(i).getNumberOfParentEdges()) {
//       rootGraphPointId = i;
//       cout << "WTF " << i << endl;
//     }
//   }
//   cout << "[ttkMandatoryCriticalPoints] rootGraphPointId = ";
//   cout << rootGraphPointId << endl;
//   vector<pair<double,double> >
// interval(mandatorySplitTree.getNumberOfPoints());
//   vector<double> xCoord(mandatorySplitTree.getNumberOfPoints());
//   vector<double> yCoord(mandatorySplitTree.getNumberOfPoints());
//   interval[rootGraphPointId].first = -2.0;
//   interval[rootGraphPointId].second = 2.0;
//   queue<int> graphPointQueue;
//   graphPointQueue.push(rootGraphPointId);
//   while(!graphPointQueue.empty()) {
//     int graphPoint = graphPointQueue.front();
//     graphPointQueue.pop();
//     int numberOfChildren =
// mandatorySplitTree.getPoint(graphPoint).getNumberOfChildEdges();
//     ////
//     xCoord[graphPoint] = (interval[graphPoint].first +
// interval[graphPoint].second) / 2.0;
//     ////
//     for(int i=0 ; i<numberOfChildren ; i++) {
//       int childGraphEdge =
// mandatorySplitTree.getPoint(graphPoint).getChildEdgeId(i);
//       int childGraphPoint =
// mandatorySplitTree.getEdge(childGraphEdge).getChildPointId();
//       cout << "  graphPoint " << graphPoint;
//       cout << "  childGraphPoint " << childGraphPoint;
//       cout << endl;
//       interval[childGraphPoint].first = interval[graphPoint].first +
// ((interval[graphPoint].second - interval[graphPoint].first) * i /
// numberOfChildren);
//       interval[childGraphPoint].second = interval[graphPoint].first +
// ((interval[graphPoint].second - interval[graphPoint].first) * (i+1) /
// numberOfChildren);
//       graphPointQueue.push(childGraphPoint);
//     }
//   }
//   // y coordinates
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfPoints() ; i++) {
//     int criticalPointId = mandatorySplitTree.getPoint(i).getId();
//     pair<double,double> criticalInterval;
//     if(mandatorySplitTree.getPoint(i).isLeaf()) {
//       criticalInterval =
// mandatoryCriticalPoints_.getMaximumInterval(criticalPointId);
//       // yCoord[i] = criticalInterval.first;
//     } else {
//       criticalInterval =
// mandatoryCriticalPoints_.getSplitSaddleInterval(criticalPointId);
//       // yCoord[i] = criticalInterval.second;
//     }
//     double y = 0.5 * (criticalInterval.first + criticalInterval.second);
//     yCoord[i] = y;
//   }
//   // x coordinates fix
//   vector<pair<int,double> >
// sortedXCoord(mandatorySplitTree.getNumberOfPoints());
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfPoints() ; i++) {
//     sortedXCoord[i].first = i;
//     sortedXCoord[i].second = xCoord[i];
//   }
//   sort(sortedXCoord.begin(), sortedXCoord.end(), xCoordCmp);
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfPoints() ; i++) {
//     xCoord[sortedXCoord[i].first] = i * (1.0 / ((int)
// mandatorySplitTree.getNumberOfPoints() - 1.0));
//   }
//
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfPoints() ; i++) {
//     cout << "  Point (" << i << ") ";
//     cout << "Interval [" << interval[i].first << ", " << interval[i].second
// << "]";
//     cout << " -> x = " << xCoord[i];
//     cout << ", y = " << yCoord[i] << endl;
//   }
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfPoints() ; i++) {
//     mandatorySplitTreePoints_->InsertNextPoint(xCoord[i], yCoord[i], 0.0);
//   }
//   outputSplitTree->SetPoints(mandatorySplitTreePoints_);
//   outputSplitTree->Allocate();
//   for(int i=0 ; i<mandatorySplitTree.getNumberOfEdges() ; i++) {
//     int parentId = mandatorySplitTree.getEdge(i).getParentPointId();
//     int childId = mandatorySplitTree.getEdge(i).getChildPointId();
//     mandatorySplitTreeEdge_[i]->InsertNextId(parentId);
//     mandatorySplitTreeEdge_[i]->InsertNextId(childId);
//     outputSplitTree->InsertNextCell(VTK_LINE, mandatorySplitTreeEdge_[i]);
//   }
// }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int ttkMandatoryCriticalPoints::buildVtkTree(
  vtkUnstructuredGrid *outputTree, MandatoryCriticalPoints::TreeType treeType) {

  /* Graph Informations */
  const Graph *graph = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
                         ? mandatoryCriticalPoints_.getJoinTreeGraph()
                         : mandatoryCriticalPoints_.getSplitTreeGraph();
  const vector<double> *xCoord
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getJoinTreeXLayout()
        : mandatoryCriticalPoints_.getSplitTreeXLayout();
  const vector<double> *yCoord
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getJoinTreeYLayout()
        : mandatoryCriticalPoints_.getSplitTreeYLayout();
  const vector<int> *mdtTreePointComponentId
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getMdtJoinTreePointComponentId()
        : mandatoryCriticalPoints_.getMdtSplitTreePointComponentId();
  const vector<MandatoryCriticalPoints::PointType> *mdtTreePointType
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getMdtJoinTreePointType()
        : mandatoryCriticalPoints_.getMdtSplitTreePointType();
  const vector<double> *mdtTreePointLowInterval
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getMdtJoinTreePointLowInterval()
        : mandatoryCriticalPoints_.getMdtSplitTreePointLowInterval();
  const vector<double> *mdtTreePointUpInterval
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getMdtJoinTreePointUpInterval()
        : mandatoryCriticalPoints_.getMdtSplitTreePointUpInterval();
  const vector<int> *mdtTreeEdgeSwitchable
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? mandatoryCriticalPoints_.getMdtJoinTreeEdgeSwitchable()
        : mandatoryCriticalPoints_.getMdtSplitTreeEdgeSwitchable();
  const int numberOfPoints = graph->getNumberOfVertices();
  const int numberOfEdges = graph->getNumberOfEdges();

  /* VTK Objects */
  vtkPoints **mdtTreePoints
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mandatoryJoinTreePoints_)
        : &(mandatorySplitTreePoints_);
  vector<vtkIdList *> *mdtTreeEdge
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mandatoryJoinTreeEdge_)
        : &(mandatorySplitTreeEdge_);
  vtkIntArray **outputTreePointType
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mdtJoinTreePointType_)
        : &(mdtSplitTreePointType_);
  vtkDoubleArray **outputTreePointLowInterval
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mdtJoinTreePointLowInterval_)
        : &(mdtSplitTreePointLowInterval_);
  vtkDoubleArray **outputTreePointUpInterval
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mdtJoinTreePointUpInterval_)
        : &(mdtSplitTreePointUpInterval_);
  vtkIntArray **outputTreePointComponentId
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mdtJoinTreePointComponentId_)
        : &(mdtSplitTreePointComponentId_);

  vtkIntArray **outputTreeEdgeSwitchable
    = (treeType == MandatoryCriticalPoints::TreeType::JoinTree)
        ? &(mdtJoinTreeEdgeSwitchable_)
        : &(mdtSplitTreeEdgeSwitchable_);

  /* Point datas */
  if((*outputTreePointType) == NULL) {
    (*outputTreePointType) = vtkIntArray::New();
    (*outputTreePointType)->SetName("Type");
  }
  (*outputTreePointType)->SetNumberOfTuples(numberOfPoints);
  outputTree->GetPointData()->AddArray(*outputTreePointType);
  if((*outputTreePointLowInterval) == NULL) {
    (*outputTreePointLowInterval) = vtkDoubleArray::New();
    (*outputTreePointLowInterval)->SetName("LowInterval");
  }
  (*outputTreePointLowInterval)->SetNumberOfTuples(numberOfPoints);
  outputTree->GetPointData()->AddArray(*outputTreePointLowInterval);
  if((*outputTreePointUpInterval) == NULL) {
    (*outputTreePointUpInterval) = vtkDoubleArray::New();
    (*outputTreePointUpInterval)->SetName("UpInterval");
  }
  (*outputTreePointUpInterval)->SetNumberOfTuples(numberOfPoints);
  outputTree->GetPointData()->AddArray(*outputTreePointUpInterval);
  if((*outputTreePointComponentId) == NULL) {
    (*outputTreePointComponentId) = vtkIntArray::New();
    (*outputTreePointComponentId)->SetName("ComponentId");
  }
  (*outputTreePointComponentId)->SetNumberOfTuples(numberOfPoints);
  outputTree->GetPointData()->AddArray(*outputTreePointComponentId);

  for(int i = 0; i < numberOfPoints; i++) {
    (*outputTreePointType)->SetTuple1(i, (*mdtTreePointType)[i]);
  }
  for(int i = 0; i < numberOfPoints; i++) {
    (*outputTreePointLowInterval)->SetTuple1(i, (*mdtTreePointLowInterval)[i]);
  }
  for(int i = 0; i < numberOfPoints; i++) {
    (*outputTreePointUpInterval)->SetTuple1(i, (*mdtTreePointUpInterval)[i]);
  }
  for(int i = 0; i < numberOfPoints; i++) {
    (*outputTreePointComponentId)->SetTuple1(i, (*mdtTreePointComponentId)[i]);
  }

  /* Cell datas */
  if((*outputTreeEdgeSwitchable) == NULL) {
    (*outputTreeEdgeSwitchable) = vtkIntArray::New();
    (*outputTreeEdgeSwitchable)->SetName("Switchable");
  }
  (*outputTreeEdgeSwitchable)->SetNumberOfTuples(numberOfEdges);
  outputTree->GetCellData()->AddArray(*outputTreeEdgeSwitchable);
  for(int i = 0; i < numberOfEdges; i++) {
    (*outputTreeEdgeSwitchable)->SetTuple1(i, (*mdtTreeEdgeSwitchable)[i]);
  }

  // Clear the VTK objects
  if((*mdtTreePoints) != NULL)
    (*mdtTreePoints)->Delete();
  if(mdtTreeEdge->size()) {
    for(size_t i = 0; i < mdtTreeEdge->size(); i++)
      if((*mdtTreeEdge)[i] != NULL)
        (*mdtTreeEdge)[i]->Delete();
    mdtTreeEdge->clear();
  }

  /* New Objects */
  // Points
  (*mdtTreePoints) = vtkPoints::New();
  outputTree->SetPoints((*mdtTreePoints));
  // Edges
  mdtTreeEdge->resize(graph->getNumberOfEdges(), NULL);
  for(size_t i = 0; i < mdtTreeEdge->size(); i++)
    (*mdtTreeEdge)[i] = vtkIdList::New();
  outputTree->Allocate(numberOfEdges);

  // Points
  for(int i = 0; i < numberOfPoints; i++) {
    (*mdtTreePoints)->InsertNextPoint((*xCoord)[i], (*yCoord)[i], 0.0);
  }

  // Edges (cells)
  for(int i = 0; i < numberOfEdges; i++) {
    (*mdtTreeEdge)[i]->InsertNextId(graph->getEdge(i)->getVertexIdx().first);
    (*mdtTreeEdge)[i]->InsertNextId(graph->getEdge(i)->getVertexIdx().second);
    outputTree->InsertNextCell(VTK_LINE, (*mdtTreeEdge)[i]);
  }
  return 0;
}

// int ttkMandatoryCriticalPoints::buildMandatoryTree(vtkUnstructuredGrid
// *outputJoinTree, const Digraph &mandatoryTree, const bool isJoinTree) {
//
//   // Clear the objects
//   if(mandatoryJoinTreePoints_)
//     mandatoryJoinTreePoints_->Delete();
//   if(mandatoryJoinTreeEdge_.size()) {
//     for(int i=0 ; i<(int)mandatoryJoinTreeEdge_.size() ; i++)
//       if(mandatoryJoinTreeEdge_[i])
//         mandatoryJoinTreeEdge_[i]->Delete();
//     mandatoryJoinTreeEdge_.clear();
//   }
//
//   if(!mandatoryTree.getNumberOfPoints())
//     return -1;
//
//   // Tree Informations
//   // const Digraph &mandatoryJoinTree =
// mandatoryCriticalPoints_.getMandatoryJoinTree();
//
//   /* New Objects */
//   // Points
//   mandatoryJoinTreePoints_ = vtkPoints::New();
//   // Edges
//   mandatoryJoinTreeEdge_.resize(mandatoryTree.getNumberOfEdges(),NULL);
//   for(int i=0 ; i<(int)mandatoryJoinTreeEdge_.size() ; i++)
//     mandatoryJoinTreeEdge_[i] = vtkIdList::New();
//
//   cout << "[ttkMandatoryCriticalPoints] Number of points : ";
//   cout << mandatoryTree.getNumberOfPoints() << endl;
//   cout << "[ttkMandatoryCriticalPoints] Number of edges : ";
//   cout << mandatoryJoinTreeEdge_.size() << endl;
//
//   // PointId to GraphPointId
//   vector<int> graphPointToVtkPoint(mandatoryTree.getNumberOfPoints(), -1);
//   // Level
//   vector<int> pointLevel(mandatoryTree.getNumberOfPoints(), -1);
//   // Find the root
//   int rootGraphPointId = -1;
//   for(int i=0 ; i<(mandatoryTree.getNumberOfPoints()) ; i++) {
//     if(!mandatoryTree.getPoint(i).getNumberOfParentEdges()) {
//       rootGraphPointId = i;
//       cout << "WTF " << i << endl;
//     }
//   }
//
//   cout << "[ttkMandatoryCriticalPoints] rootGraphPointId = ";
//   cout << rootGraphPointId << endl;
//
//   vector<pair<double,double> > interval(mandatoryTree.getNumberOfPoints());
//   vector<double> xCoord(mandatoryTree.getNumberOfPoints());
//   vector<double> yCoord(mandatoryTree.getNumberOfPoints());
//   interval[rootGraphPointId].first = -2.0;
//   interval[rootGraphPointId].second = 2.0;
//
//   queue<int> graphPointQueue;
//   graphPointQueue.push(rootGraphPointId);
//   while(!graphPointQueue.empty()) {
//     int graphPoint = graphPointQueue.front();
//     graphPointQueue.pop();
//     int numberOfChildren =
// mandatoryTree.getPoint(graphPoint).getNumberOfChildEdges();
//     ////
//     xCoord[graphPoint] = (interval[graphPoint].first +
// interval[graphPoint].second) / 2.0;
//     ////
//     for(int i=0 ; i<numberOfChildren ; i++) {
//       int childGraphEdge =
// mandatoryTree.getPoint(graphPoint).getChildEdgeId(i);
//       int childGraphPoint =
// mandatoryTree.getEdge(childGraphEdge).getChildPointId();
//       cout << "  graphPoint " << graphPoint;
//       cout << "  childGraphPoint " << childGraphPoint;
//       cout << endl;
//       interval[childGraphPoint].first = interval[graphPoint].first +
// ((interval[graphPoint].second - interval[graphPoint].first) * i /
// numberOfChildren);
//       interval[childGraphPoint].second = interval[graphPoint].first +
// ((interval[graphPoint].second - interval[graphPoint].first) * (i+1) /
// numberOfChildren);
//       graphPointQueue.push(childGraphPoint);
//     }
//   }
//   cout << "WTF Y COORDINATES" << endl;
//   // y coordinates
//   for(int i=0 ; i<mandatoryTree.getNumberOfPoints() ; i++) {
//     int criticalPointId = mandatoryTree.getPoint(i).getId();
//     pair<double,double> criticalInterval;
//     if(mandatoryTree.getPoint(i).isLeaf()) {
//       if(isJoinTree)
//         criticalInterval =
// mandatoryCriticalPoints_.getMinimumInterval(criticalPointId);
//       else
//         criticalInterval =
// mandatoryCriticalPoints_.getMaximumInterval(criticalPointId);
//       // yCoord[i] = criticalInterval.first;
//     } else {
//       if(isJoinTree)
//         criticalInterval =
// mandatoryCriticalPoints_.getJoinSaddleInterval(criticalPointId);
//       else
//         criticalInterval =
// mandatoryCriticalPoints_.getSplitSaddleInterval(criticalPointId);
//       // yCoord[i] = criticalInterval.second;
//     }
//     double y = 0.5 * (criticalInterval.first + criticalInterval.second);
//     yCoord[i] = y;
//   }
//   cout << "WTF X COORDINATES" << endl;
//   // x coordinates fix
//   if(mandatoryTree.getNumberOfPoints() > 1){
//     vector<pair<int,double> >
// sortedXCoord(mandatoryTree.getNumberOfPoints());
//     for(int i=0 ; i<mandatoryTree.getNumberOfPoints() ; i++) {
//       sortedXCoord[i].first = i;
//       sortedXCoord[i].second = xCoord[i];
//     }
//     sort(sortedXCoord.begin(), sortedXCoord.end(), xCoordCmp);
//     for(int i=0 ; i<mandatoryTree.getNumberOfPoints() ; i++) {
//       xCoord[sortedXCoord[i].first] = i * (1.0 / ((int)
// mandatoryTree.getNumberOfPoints() - 1.0));
//     }
//   }
//   cout << "WTF FIN" << endl;
//
//   for(int i=0 ; i<mandatoryTree.getNumberOfPoints() ; i++) {
//     cout << "  Point (" << i << ") ";
//     cout << "Interval [" << interval[i].first << ", " << interval[i].second
// << "]";
//     cout << " -> x = " << xCoord[i];
//     cout << ", y = " << yCoord[i] << endl;
//   }
//
//   for(int i=0 ; i<mandatoryTree.getNumberOfPoints() ; i++) {
//     mandatoryJoinTreePoints_->InsertNextPoint(xCoord[i], yCoord[i], 0.0);
//   }
//   outputJoinTree->SetPoints(mandatoryJoinTreePoints_);
//   outputJoinTree->Allocate();
//   for(int i=0 ; i<mandatoryTree.getNumberOfEdges() ; i++) {
//     int parentId = mandatoryTree.getEdge(i).getParentPointId();
//     int childId = mandatoryTree.getEdge(i).getChildPointId();
//     mandatoryJoinTreeEdge_[i]->InsertNextId(parentId);
//     mandatoryJoinTreeEdge_[i]->InsertNextId(childId);
//     outputJoinTree->InsertNextCell(VTK_LINE, mandatoryJoinTreeEdge_[i]);
//   }
//
//   return 0;
//
//
//
//
//
//   // // Points
//   // mandatoryJoinTreePoints_->InsertNextPoint(0,0,0);
//   // mandatoryJoinTreePoints_->InsertNextPoint(1,0,0);
//   // mandatoryJoinTreePoints_->InsertNextPoint(0,1,0);
//   // // Lists
//   // mandatoryJoinTreeEdge_[0]->InsertNextId(0);
//   // mandatoryJoinTreeEdge_[0]->InsertNextId(1);
//   // mandatoryJoinTreeEdge_[1]->InsertNextId(0);
//   // mandatoryJoinTreeEdge_[1]->InsertNextId(2);
//   // // Building of the structure
//   // outputJoinTree->SetPoints(mandatoryJoinTreePoints_);
//   // outputJoinTree->Allocate();
//   // outputJoinTree->InsertNextCell(VTK_LINE, mandatoryJoinTreeEdge_[0]);
//   // outputJoinTree->InsertNextCell(VTK_LINE, mandatoryJoinTreeEdge_[1]);
//   //
//   // for(int i=0 ; i<(int)minimumParentSaddleId.size() ; i++) {
//   //   if(/*minimumParentSaddleId[i] != -1*/ true) {
//   //     mandatoryJoinTreePoints_->InsertNextPoint(i/23.0,0,0);
//   //   }
//   // }
//   // outputJoinTree->SetPoints(mandatoryJoinTreePoints_);
//   // outputJoinTree->Allocate();
//   // for(int i=0 ; i<(int)minimumParentSaddleId.size() ; i++) {
//   //   mandatoryJoinTreeEdge_[i]->InsertNextId(i);
//   //   // mandatoryJoinTreeEdge_[i]->InsertNextId(i+1);
//   //   outputJoinTree->InsertNextCell(VTK_VERTEX, mandatoryJoinTreeEdge_[i]);
//   // }
// }

// transmit abort signals -- to copy paste in other wrappers
bool ttkMandatoryCriticalPoints::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkMandatoryCriticalPoints::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkMandatoryCriticalPoints] " << progress * 100
        << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkMandatoryCriticalPoints::doIt(vector<vtkDataSet *> &inputs,
                                     vector<vtkDataSet *> &outputs) {
  //   vtkDataSet *input, vtkDataSet
  // *outputMinimum, vtkDataSet *outputJoinSaddle, vtkDataSet
  // *outputSplitSaddle, vtkDataSet *outputMaximum, vtkUnstructuredGrid
  // *outputJoinTree, vtkUnstructuredGrid *outputSplitTree){

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkDataSet *outputMinimum = outputs[0];
  vtkDataSet *outputJoinSaddle = outputs[1];
  vtkDataSet *outputSplitSaddle = outputs[2];
  vtkDataSet *outputMaximum = outputs[3];
  vtkUnstructuredGrid *outputJoinTree
    = vtkUnstructuredGrid::SafeDownCast(outputs[4]);
  vtkUnstructuredGrid *outputSplitTree
    = vtkUnstructuredGrid::SafeDownCast(outputs[5]);

  // Check the last modification of the input
  if(inputMTime_ != input->GetMTime()) {
    inputMTime_ = input->GetMTime();
    computeAll_ = true;
  }

  // Use a pointer-base copy for the input data
  outputMinimum->ShallowCopy(input);
  outputJoinSaddle->ShallowCopy(input);
  outputSplitSaddle->ShallowCopy(input);
  outputMaximum->ShallowCopy(input);

  // Input data arrays
  vtkDataArray *inputLowerBoundField = NULL;
  vtkDataArray *inputUpperBoundField = NULL;

  // Get the upper bound field array in the input data set
  if(upperBoundFiledName_.length()) {
    inputUpperBoundField
      = input->GetPointData()->GetArray(upperBoundFiledName_.c_str());
  } else {
    inputUpperBoundField = input->GetPointData()->GetArray(upperBoundId);
  }
  // Error if not found
  if(!inputUpperBoundField) {
    return -1;
  }

  // Get the lower bound field array in the input data set
  if(lowerBoundFieldName_.length()) {
    inputLowerBoundField
      = input->GetPointData()->GetArray(lowerBoundFieldName_.c_str());
  } else {
    inputLowerBoundField = input->GetPointData()->GetArray(lowerBoundId);
  }
  // Error if not found
  if(!inputLowerBoundField) {
    return -1;
  }

  {
    stringstream msg;
    msg << "[ttkMandatoryCriticalPoints] Using `"
        << inputLowerBoundField->GetName() << "' as lower bound..." << endl;
    msg << "[ttkMandatoryCriticalPoints] Using `"
        << inputUpperBoundField->GetName() << "' as upper bound..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Initialize the triangulation object with the input data set
  triangulation_ = ttkTriangulation::getTriangulation(input);

  if(!triangulation_)
    return -1;

  triangulation_->setWrapper(this);
  mandatoryCriticalPoints_.setupTriangulation(triangulation_);

  bool hasChangedConnectivity = false;

  if((triangulation_->isEmpty())
     || (ttkTriangulation::hasChangedConnectivity(
          triangulation_, input, this))) {
    Modified();
    hasChangedConnectivity = true;
  }

  // Allocate the memory for the output scalar field
  if((!outputMandatoryMinimum_) || (hasChangedConnectivity)) {
    if(!outputMandatoryMinimum_)
      outputMandatoryMinimum_ = vtkIntArray::New();
    outputMandatoryMinimum_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputMandatoryMinimum_->SetName("MinimumComponents");
  }
  if((!outputMandatoryJoinSaddle_) || (hasChangedConnectivity)) {
    if(!outputMandatoryJoinSaddle_)
      outputMandatoryJoinSaddle_ = vtkIntArray::New();
    outputMandatoryJoinSaddle_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputMandatoryJoinSaddle_->SetName("JoinSaddleComponents");
  }
  if((!outputMandatorySplitSaddle_) || (hasChangedConnectivity)) {
    if(!outputMandatorySplitSaddle_)
      outputMandatorySplitSaddle_ = vtkIntArray::New();
    outputMandatorySplitSaddle_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputMandatorySplitSaddle_->SetName("SplitSaddleComponents");
  }
  if((!outputMandatoryMaximum_) || (hasChangedConnectivity)) {
    if(!outputMandatoryMaximum_)
      outputMandatoryMaximum_ = vtkIntArray::New();
    outputMandatoryMaximum_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputMandatoryMaximum_->SetName("MaximumComponents");
  }
  // Add array in the output data set
  outputMinimum->GetPointData()->AddArray(outputMandatoryMinimum_);
  outputJoinSaddle->GetPointData()->AddArray(outputMandatoryJoinSaddle_);
  outputSplitSaddle->GetPointData()->AddArray(outputMandatorySplitSaddle_);
  outputMaximum->GetPointData()->AddArray(outputMandatoryMaximum_);

  // Calling the executing package
  switch(inputUpperBoundField->GetDataType()) {

    vtkTemplateMacro({
      // Reset the baseCode object
      if((computeAll_) || (hasChangedConnectivity))
        mandatoryCriticalPoints_.flush();

      // Wrapper
      mandatoryCriticalPoints_.setWrapper(this);
      // Set the number of vertex
      mandatoryCriticalPoints_.setVertexNumber(input->GetNumberOfPoints());
      // Set the coordinates of each vertex
      double point[3];
      for(int i = 0; i < input->GetNumberOfPoints(); i++) {
        input->GetPoint(i, point);
        mandatoryCriticalPoints_.setVertexPosition(i, point);
      }
      // Set the void pointers to the upper and lower bound fields
      mandatoryCriticalPoints_.setLowerBoundFieldPointer(
        inputLowerBoundField->GetVoidPointer(0));
      mandatoryCriticalPoints_.setUpperBoundFieldPointer(
        inputUpperBoundField->GetVoidPointer(0));
      // Set the output data pointers
      mandatoryCriticalPoints_.setOutputMinimumDataPointer(
        outputMandatoryMinimum_->GetVoidPointer(0));
      mandatoryCriticalPoints_.setOutputJoinSaddleDataPointer(
        outputMandatoryJoinSaddle_->GetVoidPointer(0));
      mandatoryCriticalPoints_.setOutputSplitSaddleDataPointer(
        outputMandatorySplitSaddle_->GetVoidPointer(0));
      mandatoryCriticalPoints_.setOutputMaximumDataPointer(
        outputMandatoryMaximum_->GetVoidPointer(0));
      // Set the triangulation object
      // Set the offsets
      mandatoryCriticalPoints_.setSoSoffsets();
      // Simplification threshold

      mandatoryCriticalPoints_.setSimplificationThreshold(
        simplificationThreshold_);
      // Execute process
      if(computeAll_) {
        mandatoryCriticalPoints_.execute<VTK_TT>();
        computeAll_ = false;
        simplify_ = false;
        computeMinimumOutput_ = true;
        computeJoinSaddleOutput_ = true;
        computeSplitSaddleOutput_ = true;
        computeMaximumOutput_ = true;
      }
    });
  }

  // Simplification
  if(simplify_) {
    // Simplify
    mandatoryCriticalPoints_.simplifyJoinTree();
    mandatoryCriticalPoints_.buildJoinTreePlanarLayout();
    mandatoryCriticalPoints_.simplifySplitTree();
    mandatoryCriticalPoints_.buildSplitTreePlanarLayout();
    // Recompute the outputs
    computeMinimumOutput_ = true;
    computeJoinSaddleOutput_ = true;
    computeSplitSaddleOutput_ = true;
    computeMaximumOutput_ = true;
  }
  // Outputs
  if(computeMinimumOutput_) {
    if(outputAllMinimumComponents_)
      mandatoryCriticalPoints_.outputAllMinima();
    else
      mandatoryCriticalPoints_.outputMinimum(outputMinimumComponentId_);
    computeMinimumOutput_ = false;
  }
  if(computeJoinSaddleOutput_) {
    if(outputAllJoinSaddleComponents_)
      mandatoryCriticalPoints_.outputAllJoinSaddle();
    else
      mandatoryCriticalPoints_.outputJoinSaddle(outputJoinSaddleComponentId_);
    computeJoinSaddleOutput_ = false;
  }
  if(computeSplitSaddleOutput_) {
    if(outputAllSplitSaddleComponents_)
      mandatoryCriticalPoints_.outputAllSplitSaddle();
    else
      mandatoryCriticalPoints_.outputSplitSaddle(outputSplitSaddleComponentId_);
    computeSplitSaddleOutput_ = false;
  }
  if(computeMaximumOutput_) {
    if(outputAllMaximumComponents_)
      mandatoryCriticalPoints_.outputAllMaxima();
    else
      mandatoryCriticalPoints_.outputMaximum(outputMaximumComponentId_);
    computeMaximumOutput_ = false;
  }

  buildVtkTree(outputJoinTree, MandatoryCriticalPoints::TreeType::JoinTree);
  buildVtkTree(outputSplitTree, MandatoryCriticalPoints::TreeType::SplitTree);

  {
    stringstream msg;
    memoryUsage_ += m.getElapsedUsage();
    msg << "[ttkMandatoryCriticalPoints] Memory usage: " << memoryUsage_
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkMandatoryCriticalPoints::FillInputPortInformation(int port,
                                                         vtkInformation *info) {
  if(!this->Superclass::FillInputPortInformation(port, info))
    return 0;
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkMandatoryCriticalPoints::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(!this->Superclass::FillOutputPortInformation(port, info))
    return 0;
  if(port >= 0 && port < 4)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if((port == 4) || (port == 5))
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}
