#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include "ttkContourForests.h"

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;
using namespace cf;

vtkStandardNewMacro(ttkContourForests);

ttkContourForests::ttkContourForests() {
  // VTK Interface //
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);

  vtkWarningMacro(
    "Contour Forests is deprecated, please use FTM Tree instead.");
}

void ttkContourForests::clearSkeleton() {
  samples_.clear();
  barycenters_.clear();

  skeletonNodes_->Delete();
  skeletonNodes_ = vtkPolyData::New();
  skeletonArcs_->Delete();
  skeletonArcs_ = vtkPolyData::New();
}

void ttkContourForests::clearSegmentation() {
  segmentation_->Delete();
}

void ttkContourForests::clearTree() {
  tree_ = nullptr;
}

int ttkContourForests::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkContourForests::FillOutputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void ttkContourForests::Modified() {
  toComputeSkeleton_ = true;
  toComputeContourTree_ = true;
  ttkAlgorithm::Modified();
}

void ttkContourForests::SetForceInputOffsetScalarField(bool onOff) {
  toUpdateVertexSoSoffsets_ = true;
  toComputeContourTree_ = true;
  toUpdateTree_ = true;
  toComputeSkeleton_ = true;
  toComputeSegmentation_ = true;

  ForceInputOffsetScalarField = onOff;
  Modified();
}

void ttkContourForests::SetTreeType(int treeType) {
  if(treeType >= 0 && treeType <= 2) {
    toUpdateTree_ = true;
    if(treeType_ == TreeType::Contour
       || static_cast<TreeType>(treeType) == TreeType::Contour) {
      toComputeContourTree_ = true;
    }

    toComputeSkeleton_ = true;
    toComputeSegmentation_ = true;

    treeType_ = static_cast<TreeType>(treeType);

    Modified();
  }
}

void ttkContourForests::ShowMin(bool state) {
  toComputeSkeleton_ = true;

  showMin_ = state;
  Modified();
}

void ttkContourForests::ShowMax(bool state) {
  toComputeSkeleton_ = true;

  showMax_ = state;
  Modified();
}

void ttkContourForests::ShowSaddle1(bool state) {
  toComputeSkeleton_ = true;

  showSaddle1_ = state;
  Modified();
}

void ttkContourForests::ShowSaddle2(bool state) {
  toComputeSkeleton_ = true;

  showSaddle2_ = state;
  Modified();
}

void ttkContourForests::ShowArc(bool state) {
  toComputeSkeleton_ = true;

  showArc_ = state;
  Modified();
}

void ttkContourForests::SetArcResolution(int arcResolution) {
  if(arcResolution >= 0) {
    toComputeSkeleton_ = true;

    arcResolution_ = arcResolution;
    Modified();
  }
}

void ttkContourForests::SetPartitionNumber(int partitionNum) {
  partitionNum_ = partitionNum;

  toComputeContourTree_ = true;
  toComputeSkeleton_ = true;
  toUpdateTree_ = true;
  Modified();
}

void ttkContourForests::SetLessPartition(bool l) {
  lessPartition_ = l;
  Modified();
}

void ttkContourForests::SetSkeletonSmoothing(double skeletonSmoothing) {
  if(skeletonSmoothing >= 0) {
    toComputeSkeleton_ = true;

    skeletonSmoothing_ = skeletonSmoothing;
    Modified();
  }
}

void ttkContourForests::SetSimplificationType(int type) {
  simplificationType_ = type;
  Modified();
}

void ttkContourForests::SetSimplificationThreshold(
  double simplificationThreshold) {
  if(simplificationThreshold >= 0.0 && simplificationThreshold <= 1.0) {
    simplificationThresholdBuffer_ = simplificationThreshold;
    toComputeContourTree_ = true;
    toUpdateTree_ = true;
    toComputeSkeleton_ = true;
    toComputeSegmentation_ = true;

    Modified();
  }
}

bool ttkContourForests::isCoincident(float p1[], double p2[]) {
  double sPrev[3];
  double sNext[3];
  for(unsigned int k = 0; k < 3; k++) {
    sPrev[k] = p2[k] - p1[k];
    sNext[k] = sPrev[k];
  }

  return (vtkMath::Normalize(sNext) == 0.0);
}

bool ttkContourForests::isCoincident(double p1[], double p2[]) {
  double sPrev[3];
  double sNext[3];
  for(unsigned int k = 0; k < 3; k++) {
    sPrev[k] = p2[k] - p1[k];
    sNext[k] = sPrev[k];
  }

  return (vtkMath::Normalize(sNext) == 0.0);
}

void ttkContourForests::getSkeletonArcs() {
  vtkNew<vtkAppendPolyData> app{};

  vtkDoubleArray *scalars{};
  ttkSimplexIdTypeArray *identifierScalars{};
  vtkIntArray *typeScalars{};
  ttkSimplexIdTypeArray *sizeScalars{};
  vtkDoubleArray *spanScalars{};
  int type = static_cast<int>(TreeComponent::Arc);

  float point1[3];
  vector<double> point2(3);
  // get skeleton scalars
  vector<vector<double>> skeletonScalars{};
  getSkeletonScalars(vertexScalars_, skeletonScalars);

  double inputScalar;
  SuperArc *a;
  SimplexId regionSize;
  double regionSpan;
  SimplexId currentZone = 0;
  SimplexId regionId;

  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if(a->isVisible()) {
      SimplexId upNodeId = tree_->getSuperArc(i)->getUpNodeId();
      SimplexId upVertex = tree_->getNode(upNodeId)->getVertexId();
      float coordUp[3];
      triangulation_->getVertexPoint(
        upVertex, coordUp[0], coordUp[1], coordUp[2]);

      SimplexId downNodeId = tree_->getSuperArc(i)->getDownNodeId();
      SimplexId downVertex = tree_->getNode(downNodeId)->getVertexId();
      float coordDown[3];
      triangulation_->getVertexPoint(
        downVertex, coordDown[0], coordDown[1], coordDown[2]);

      regionSize = tree_->getNumberOfVisibleRegularNode(i);
      regionSpan = Geometry::distance(coordUp, coordDown, 3);
      regionId = currentZone++;

      // Line //
      if(barycenters_[static_cast<int>(treeType_)][i].size()) {
        // init: min
        SimplexId downNodeVId;
        if(treeType_ == TreeType::Split)
          downNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
        else
          downNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();

        triangulation_->getVertexPoint(
          downNodeVId, point1[0], point1[1], point1[2]);
        vtkSmartPointer<vtkLineSource> line
          = vtkSmartPointer<vtkLineSource>::New();
        line->SetPoint1(point1);

        const auto nbBarycenter
          = barycenters_[static_cast<int>(treeType_)][i].size();
        for(unsigned int j = 0; j < nbBarycenter; ++j) {
          point2 = barycenters_[static_cast<int>(treeType_)][i][j];
          line->SetPoint2(point2.data());

          if(!isCoincident(point1, point2.data())) {
            line->Update();
            vtkPolyData *lineData = line->GetOutput();

            // Point data //
            {
              inputScalar = skeletonScalars[i][j];

              scalars = vtkDoubleArray::New();
              scalars->SetName(vtkInputScalars_->GetName());
              for(unsigned int k = 0; k < 2; ++k)
                scalars->InsertTuple1(k, inputScalar);
              lineData->GetPointData()->AddArray(scalars);
              scalars->Delete();
            }

            // Cell data //
            // Identifier
            identifierScalars = ttkSimplexIdTypeArray::New();
            identifierScalars->SetName("SegmentationId");
            for(unsigned int k = 0; k < 2; ++k)
              identifierScalars->InsertTuple1(k, regionId);
            lineData->GetCellData()->AddArray(identifierScalars);
            identifierScalars->Delete();
            // Type
            typeScalars = vtkIntArray::New();
            typeScalars->SetName("Type");
            for(unsigned int k = 0; k < 2; ++k)
              typeScalars->InsertTuple1(k, type);
            lineData->GetCellData()->AddArray(typeScalars);
            typeScalars->Delete();
            // Size
            sizeScalars = ttkSimplexIdTypeArray::New();
            sizeScalars->SetName("RegionSize");
            for(unsigned int k = 0; k < 2; ++k)
              sizeScalars->InsertTuple1(k, regionSize);
            lineData->GetCellData()->AddArray(sizeScalars);
            sizeScalars->Delete();
            // Span
            spanScalars = vtkDoubleArray::New();
            spanScalars->SetName("RegionSpan");
            for(unsigned int k = 0; k < 2; ++k)
              spanScalars->InsertTuple1(k, regionSpan);
            lineData->GetCellData()->AddArray(spanScalars);
            spanScalars->Delete();

            app->AddInputData(lineData);
          }

          line = vtkSmartPointer<vtkLineSource>::New();
          line->SetPoint1(point2.data());
        }

        // end: max
        SimplexId upNodeVId;
        if(treeType_ == TreeType::Split)
          upNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();
        else
          upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();

        std::array<float, 3> pt{};
        triangulation_->getVertexPoint(upNodeVId, pt[0], pt[1], pt[2]);
        point2[0] = pt[0];
        point2[1] = pt[1];
        point2[2] = pt[2];
        line->SetPoint2(point2.data());

        if(!isCoincident(point1, point2.data())) {
          line->Update();
          vtkPolyData *lineData = line->GetOutput();

          // Point data //
          {
            inputScalar
              = skeletonScalars[i][barycenters_[static_cast<int>(treeType_)][i]
                                     .size()];

            scalars = vtkDoubleArray::New();
            scalars->SetName(vtkInputScalars_->GetName());
            for(unsigned int k = 0; k < 2; ++k)
              scalars->InsertTuple1(k, inputScalar);
            lineData->GetPointData()->AddArray(scalars);
            scalars->Delete();
          }

          // Cell data //
          // Identifier
          identifierScalars = ttkSimplexIdTypeArray::New();
          identifierScalars->SetName("SegmentationId");
          for(unsigned int k = 0; k < 2; ++k)
            identifierScalars->InsertTuple1(k, regionId);
          lineData->GetCellData()->AddArray(identifierScalars);
          identifierScalars->Delete();
          // Type
          typeScalars = vtkIntArray::New();
          typeScalars->SetName("Type");
          for(unsigned int k = 0; k < 2; ++k)
            typeScalars->InsertTuple1(k, type);
          lineData->GetCellData()->AddArray(typeScalars);
          typeScalars->Delete();
          // Size
          sizeScalars = ttkSimplexIdTypeArray::New();
          sizeScalars->SetName("RegionSize");
          for(unsigned int k = 0; k < 2; ++k)
            sizeScalars->InsertTuple1(k, regionSize);
          lineData->GetCellData()->AddArray(sizeScalars);
          sizeScalars->Delete();
          // Span
          spanScalars = vtkDoubleArray::New();
          spanScalars->SetName("RegionSpan");
          for(unsigned int k = 0; k < 2; ++k)
            spanScalars->InsertTuple1(k, regionSpan);
          lineData->GetCellData()->AddArray(spanScalars);
          spanScalars->Delete();

          app->AddInputData(lineData);
        }
      } else {
        vtkNew<vtkLineSource> line{};

        SimplexId downNodeVId
          = tree_->getNode(a->getDownNodeId())->getVertexId();
        triangulation_->getVertexPoint(
          downNodeVId, point1[0], point1[1], point1[2]);
        line->SetPoint1(point1);

        SimplexId upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
        std::array<float, 3> pt{};
        triangulation_->getVertexPoint(upNodeVId, pt[0], pt[1], pt[2]);
        point2[0] = pt[0];
        point2[1] = pt[1];
        point2[2] = pt[2];
        line->SetPoint2(point2.data());

        if(!isCoincident(point1, point2.data())) {
          line->Update();
          vtkPolyData *lineData = line->GetOutput();

          // Point data //
          {
            inputScalar = skeletonScalars[i][0];

            scalars = vtkDoubleArray::New();
            scalars->SetName(vtkInputScalars_->GetName());
            for(unsigned int k = 0; k < 2; ++k)
              scalars->InsertTuple1(k, inputScalar);
            lineData->GetPointData()->AddArray(scalars);
            scalars->Delete();
          }

          // Cell data //
          // Identifier
          identifierScalars = ttkSimplexIdTypeArray::New();
          identifierScalars->SetName("SegmentationId");
          for(int k = 0; k < 2; ++k)
            identifierScalars->InsertTuple1(k, regionId);
          lineData->GetCellData()->AddArray(identifierScalars);
          identifierScalars->Delete();
          // Type
          typeScalars = vtkIntArray::New();
          typeScalars->SetName("Type");
          for(unsigned int k = 0; k < 2; ++k)
            typeScalars->InsertTuple1(k, type);
          lineData->GetCellData()->AddArray(typeScalars);
          typeScalars->Delete();
          // Size
          sizeScalars = ttkSimplexIdTypeArray::New();
          sizeScalars->SetName("RegionSize");
          for(unsigned int k = 0; k < 2; ++k)
            sizeScalars->InsertTuple1(k, regionSize);
          lineData->GetCellData()->AddArray(sizeScalars);
          sizeScalars->Delete();
          // Span
          spanScalars = vtkDoubleArray::New();
          spanScalars->SetName("RegionSpan");
          for(unsigned int k = 0; k < 2; ++k)
            spanScalars->InsertTuple1(k, regionSpan);
          lineData->GetCellData()->AddArray(spanScalars);
          spanScalars->Delete();

          app->AddInputData(lineData);
        }
      }
    } else {
      // cout << " pruned _ :" <<
      // tree_->getNode(a->getDownNodeId())->getVertexId() << " -  "
      //<< tree_->getNode(a->getUpNodeId())->getVertexId() << endl;
    }
  }

  app->Update();
  skeletonArcs_->ShallowCopy(app->GetOutput());
}

int ttkContourForests::getSkeletonScalars(
  const vector<double> &scalars,
  vector<vector<double>> &skeletonScalars) const {
  skeletonScalars.clear();
  skeletonScalars.resize(tree_->getNumberOfSuperArcs());

  SimplexId nodeId;
  SimplexId vertexId;

  double f;
  double f0;
  double f1;
  double fmin;
  double fmax;
  SimplexId nodeMinId;
  SimplexId nodeMaxId;
  SimplexId nodeMinVId;
  SimplexId nodeMaxVId;
  const SuperArc *a;
  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if(!a->isPruned()) {
      if(treeType_ == TreeType::Split) {
        nodeMinId = a->getUpNodeId();
        nodeMaxId = a->getDownNodeId();
      } else {
        nodeMaxId = a->getUpNodeId();
        nodeMinId = a->getDownNodeId();
      }

      nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
      nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

      fmax = scalars[nodeMaxVId];
      fmin = scalars[nodeMinVId];

      // init: min
      f0 = fmin;

      // iteration
      for(SimplexId j = 0;
          j < (SimplexId)samples_[static_cast<int>(treeType_)][i].size(); ++j) {
        const vector<SimplexId> &sample
          = samples_[static_cast<int>(treeType_)][i][j];

        f = 0;
        for(SimplexId k = 0; k < (SimplexId)sample.size(); ++k) {
          nodeId = sample[k];
          vertexId = nodeId;
          f += scalars[vertexId];
        }
        if(sample.size()) {
          f /= sample.size();

          f1 = f;
          // update the arc
          skeletonScalars[i].push_back((f0 + f1) / 2);
          f0 = f1;
        }
      }

      // end: max
      f1 = fmax;

      // update the arc
      skeletonScalars[i].push_back((f0 + f1) / 2);
    }
  }

  return 0;
}

void ttkContourForests::getSkeletonNodes() {
  vtkNew<vtkPoints> points{};
  float point[3];

  double scalar{};
  vtkNew<vtkDoubleArray> scalars{};
  scalars->SetName(vtkInputScalars_->GetName());

  vtkNew<ttkSimplexIdTypeArray> nodeIdentifierScalars{};
  nodeIdentifierScalars->SetName("NodeIdentifier");

  vtkNew<ttkSimplexIdTypeArray> vertexIdentifierScalars{};
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  int type{};
  vtkNew<vtkIntArray> nodeTypeScalars{};
  nodeTypeScalars->SetName("CriticalType");

  vtkNew<ttkSimplexIdTypeArray> regionSizeScalars{};
  regionSizeScalars->SetName("RegionSize");

  SimplexId identifier{};
  for(unsigned i = 0; i < criticalPoints_.size(); ++i) {
    SimplexId nodeId = criticalPoints_[i];
    if(tree_->getNode(nodeId)->isHidden())
      continue;
    SimplexId vertexId = tree_->getNode(nodeId)->getVertexId();
    CriticalType nodeType = getNodeType(nodeId);

    if((nodeType == CriticalType::Local_minimum and showMin_)
       or (nodeType == CriticalType::Local_maximum and showMax_)
       or (nodeType == CriticalType::Saddle1 and showSaddle1_)
       or (nodeType == CriticalType::Saddle2 and showSaddle2_)
       or ((showSaddle1_ and showSaddle2_)
           and (nodeType == CriticalType::Regular
                or nodeType == CriticalType::Degenerate))) {
      // Positions
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);
      points->InsertPoint(identifier, point);

      // Scalars
      scalar = vertexScalars_[vertexId];
      scalars->InsertTuple1(identifier, scalar);

      // NodeIdentifier
      nodeIdentifierScalars->InsertTuple1(identifier, nodeId);

      // VertexIdentifier
      vertexIdentifierScalars->InsertTuple1(identifier, vertexId);

      // Type
      type = static_cast<int>(nodeType);
      nodeTypeScalars->InsertTuple1(identifier, type);

      // RegionSize
      SimplexId regionSize = 0;
      if(nodeType == CriticalType::Local_maximum) {
        const SimplexId arcId = tree_->getNode(nodeId)->getDownSuperArcId(0);
        regionSize = tree_->getSuperArc(arcId)->getNumberOfRegularNodes() + 1;
      } else if(nodeType == CriticalType::Local_minimum) {
        const SimplexId arcId = tree_->getNode(nodeId)->getUpSuperArcId(0);
        regionSize = tree_->getSuperArc(arcId)->getNumberOfRegularNodes() + 1;
      }
      regionSizeScalars->InsertTuple1(identifier, regionSize);

      ++identifier;
    }
  }
  skeletonNodes_->SetPoints(points);
  skeletonNodes_->GetPointData()->AddArray(scalars);
  skeletonNodes_->GetPointData()->AddArray(nodeIdentifierScalars);
  skeletonNodes_->GetPointData()->AddArray(vertexIdentifierScalars);
  skeletonNodes_->GetPointData()->AddArray(nodeTypeScalars);
  skeletonNodes_->GetPointData()->AddArray(regionSizeScalars);
}

CriticalType ttkContourForests::getNodeType(SimplexId id) {
  return getNodeType(id, treeType_, tree_);
}

CriticalType
  ttkContourForests::getNodeType(SimplexId id, TreeType type, MergeTree *tree) {
  int upDegree{};
  int downDegree{};
  if(type == TreeType::Join || type == TreeType::Contour) {
    upDegree = tree->getNode(id)->getUpValence();
    downDegree = tree->getNode(id)->getDownValence();
  } else {
    downDegree = tree->getNode(id)->getUpValence();
    upDegree = tree->getNode(id)->getDownValence();
  }
  int degree = upDegree + downDegree;

  // saddle point
  if(degree > 1) {
    if(upDegree == 2 && downDegree == 1)
      return CriticalType::Saddle2;
    else if(upDegree == 1 && downDegree == 2)
      return CriticalType::Saddle1;
    else if(upDegree == 1 && downDegree == 1)
      return CriticalType::Regular;
    else
      return CriticalType::Degenerate;
  }
  // local extremum
  else {
    if(upDegree)
      return CriticalType::Local_minimum;
    else
      return CriticalType::Local_maximum;
  }
}

void ttkContourForests::getCriticalPoints() {
  vector<bool> isCriticalPoint(numberOfVertices_);

  criticalPoints_.clear();
  for(SimplexId i = 0; i < numberOfVertices_; ++i)
    isCriticalPoint[i] = false;

  // const int nbVert = triangulation_->getNumberOfOriginalVertices();

  // looking for critical points
  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    auto a = tree_->getSuperArc(i);

    if(!a->isPruned()) {
      SimplexId upId = a->getUpNodeId();
      SimplexId up_vId = tree_->getNode(upId)->getVertexId();
      if(!isCriticalPoint[up_vId]) {
        isCriticalPoint[up_vId] = true;
        criticalPoints_.push_back(upId);
      }

      SimplexId downId = a->getDownNodeId();
      SimplexId down_vId = tree_->getNode(downId)->getVertexId();
      if(!isCriticalPoint[down_vId]) {
        isCriticalPoint[down_vId] = true;
        criticalPoints_.push_back(downId);
      }
    }
  }
  //{
  // stringstream msg;
  // msg << "[ttkContourForests] List of critical points :" << endl;
  // for (unsigned int it = 0; it < criticalPoints_->size(); ++it)
  // msg << "[ttkContourForests]   NodeId:" << (*criticalPoints_)[it]
  //<< ", VertexId:" << tree_->getNode(it)->getVertexId() << endl;
  // dMsg(cout, msg.str(), advancedInfoMsg);
  //}
}

int ttkContourForests::sample(unsigned int samplingLevel) {
  samples_.resize(3);
  samples_[static_cast<int>(treeType_)].resize(tree_->getNumberOfSuperArcs());
  vector<vector<SimplexId>> sampleList(samplingLevel);

  SuperArc *a;
  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if(!a->isPruned()) {
      for(unsigned int j = 0; j < samplingLevel; ++j) {
        sampleList[j].clear();
      }

      double fmax, fmin;
      SimplexId nodeMaxId, nodeMinId;
      SimplexId nodeMaxVId, nodeMinVId;
      double delta;
      if(a->getNumberOfRegularNodes()) {
        if(treeType_ == TreeType::Split) {
          nodeMaxId = a->getDownNodeId();
          nodeMinId = a->getUpNodeId();
        } else {
          nodeMaxId = a->getUpNodeId();
          nodeMinId = a->getDownNodeId();
        }

        nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
        nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

        fmax = vertexScalars_[nodeMaxVId];
        fmin = vertexScalars_[nodeMinVId];

        delta = (fmax - fmin) / samplingLevel;

        double f;
        SimplexId nodeId;
        SimplexId vertexId;
        for(SimplexId j = 0; j < a->getNumberOfRegularNodes(); ++j) {
          nodeId = a->getRegularNodeId(j);
          if(a->isMasqued(j))
            continue;
          vertexId = nodeId;
          f = vertexScalars_[vertexId];

          for(unsigned int k = 0; k < samplingLevel; ++k) {
            if(f <= (k + 1) * delta + fmin) {
              sampleList[k].push_back(nodeId);
              break;
            }
          }
        }

        // update the arc
        for(SimplexId j = 0; j < (SimplexId)sampleList.size(); ++j)
          samples_[static_cast<int>(treeType_)][i].push_back(sampleList[j]);
      }
    }
  }

  return 0;
}

int ttkContourForests::computeBarycenters() {
  barycenters_.resize(3);
  barycenters_[static_cast<int>(treeType_)].resize(
    tree_->getNumberOfSuperArcs());
  vector<float> barycenter(3);
  SimplexId vertexId;

  const SuperArc *a;
  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);
    if(!a->isPruned()) {
      for(SimplexId j = 0;
          j < (SimplexId)samples_[static_cast<int>(treeType_)][i].size(); ++j) {
        vector<SimplexId> &sample = samples_[static_cast<int>(treeType_)][i][j];

        for(unsigned int k = 0; k < 3; ++k)
          barycenter[k] = 0;

        for(SimplexId k = 0; k < (SimplexId)sample.size(); ++k) {
          vertexId = sample[k];

          float pt[3];
          triangulation_->getVertexPoint(vertexId, pt[0], pt[1], pt[2]);
          barycenter[0] += pt[0];
          barycenter[1] += pt[1];
          barycenter[2] += pt[2];
        }
        if(sample.size()) {
          for(unsigned int k = 0; k < 3; ++k)
            barycenter[k] /= sample.size();

          // update the arc
          unsigned int nbBar
            = barycenters_[static_cast<int>(treeType_)][i].size();
          barycenters_[static_cast<int>(treeType_)][i].resize(nbBar + 1);
          barycenters_[static_cast<int>(treeType_)][i][nbBar].resize(3);

          for(unsigned int k = 0; k < 3; ++k)
            barycenters_[static_cast<int>(treeType_)][i][nbBar][k]
              = barycenter[k];
        }
      }
    }
  }

  return 0;
}

void ttkContourForests::computeSkeleton(unsigned int arcRes) {
  sample(arcRes);
  computeBarycenters();
}

void ttkContourForests::smoothSkeleton(unsigned int skeletonSmoothing) {
  for(unsigned int i = 0; i < skeletonSmoothing; i++) {
    for(SimplexId j = 0; j < (SimplexId)tree_->getNumberOfSuperArcs(); j++) {
      if(!tree_->getSuperArc(j)->isPruned()) {
        smooth(j, !(treeType_ == TreeType::Split));
      }
    }
  }
}

void ttkContourForests::smooth(const SimplexId idArc, bool order) {
  int N = barycenters_[static_cast<int>(treeType_)][idArc].size();
  if(N) {
    // init //
    vector<vector<double>> barycenterList(N);
    for(unsigned int i = 0; i < barycenterList.size(); ++i)
      barycenterList[i].resize(3);

    SimplexId up_vId;
    SimplexId down_vId;
    if(order) {
      up_vId = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())
                 ->getVertexId();
      down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())
                   ->getVertexId();
    } else {
      down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())
                   ->getVertexId();
      up_vId = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())
                 ->getVertexId();
    }

    std::array<float, 3> p0{};
    std::array<float, 3> p1{};
    triangulation_->getVertexPoint(down_vId, p0[0], p0[1], p0[2]);
    triangulation_->getVertexPoint(up_vId, p1[0], p1[1], p1[2]);

    // filtering //
    if(N > 1) {
      // first
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k]
          = (p0[k] + barycenters_[static_cast<int>(treeType_)][idArc][1][k])
            * 0.5;

      // main
      for(int i = 1; i < N - 1; ++i) {
        for(unsigned int k = 0; k < 3; ++k)
          barycenterList[i][k]
            = (barycenters_[static_cast<int>(treeType_)][idArc][i - 1][k]
               + barycenters_[static_cast<int>(treeType_)][idArc][i + 1][k])
              * 0.5;
      }
      // last
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[N - 1][k]
          = (p1[k] + barycenters_[static_cast<int>(treeType_)][idArc][N - 1][k])
            * 0.5;
    } else {
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k] = (p0[k] + p1[k]) * 0.5;
    }

    // copy //
    for(int i = 0; i < N; ++i) {
      for(unsigned int k = 0; k < 3; ++k)
        barycenters_[static_cast<int>(treeType_)][idArc][i][k]
          = barycenterList[i][k];
    }
  }
}

void ttkContourForests::getSkeleton() {
  Timer t;
  computeSkeleton(arcResolution_);
  smoothSkeleton(skeletonSmoothing_);

  // nodes
  if(showMin_ || showMax_ || showSaddle1_ || showSaddle2_)
    getSkeletonNodes();
  else
    skeletonNodes_->ShallowCopy(voidUnstructuredGrid_);

  // arcs
  if(showArc_)
    getSkeletonArcs();
  else
    skeletonArcs_->ShallowCopy(voidPolyData_);

  // what is done is no longer to be done
  toComputeSkeleton_ = false;

  this->printMsg(
    "Topological skeleton built", 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg(std::vector<std::vector<std::string>>{
    {"Arc resolution", std::to_string(arcResolution_)},
    {"Smoothing", std::to_string(skeletonSmoothing_)}});
}

void ttkContourForests::getSegmentation(vtkDataSet *input) {
  Timer t;

  // field
  SimplexId regionId{};
  vtkNew<ttkSimplexIdTypeArray> scalarsRegionId{};
  scalarsRegionId->SetName("SegmentationId");
  scalarsRegionId->SetNumberOfTuples(vertexScalars_.size());

  int regionType{};
  vtkNew<vtkIntArray> scalarsRegionType{};
  scalarsRegionType->SetName("RegionType");
  scalarsRegionType->SetNumberOfTuples(vertexScalars_.size());

  SimplexId regionSize{};
  vtkNew<ttkSimplexIdTypeArray> scalarsRegionSize{};
  scalarsRegionSize->SetName("RegionSize");
  scalarsRegionSize->SetNumberOfTuples(vertexScalars_.size());

  double regionSpan{};
  vtkNew<vtkDoubleArray> scalarsRegionSpan{};
  scalarsRegionSpan->SetName("RegionSpan");
  scalarsRegionSpan->SetNumberOfTuples(vertexScalars_.size());

  SimplexId currentZone{};

  if(!segmentation_) {
    segmentation_ = input->NewInstance();
    segmentation_->ShallowCopy(input);
  }

  for(SimplexId i = 0; i < numberOfVertices_; i++) {
    scalarsRegionId->SetTuple1(i, -1);
  }

  // nodes
  for(SimplexId it = 0; it < (SimplexId)criticalPoints_.size(); ++it) {
    SimplexId nodeId = criticalPoints_[it];
    SimplexId vertexId = tree_->getNode(nodeId)->getVertexId();

    // RegionType
    regionType = -1;
    scalarsRegionType->SetTuple1(vertexId, regionType);
  }

  // arcs
  for(SimplexId i = 0; i < (SimplexId)tree_->getNumberOfSuperArcs(); ++i) {
    auto a = tree_->getSuperArc(i);
    if(a->isVisible()) {
      SimplexId upNodeId = tree_->getSuperArc(i)->getUpNodeId();
      CriticalType upNodeType = getNodeType(upNodeId);
      SimplexId upVertex = tree_->getNode(upNodeId)->getVertexId();
      float coordUp[3];
      triangulation_->getVertexPoint(
        upVertex, coordUp[0], coordUp[1], coordUp[2]);

      SimplexId downNodeId = tree_->getSuperArc(i)->getDownNodeId();
      CriticalType downNodeType = getNodeType(downNodeId);
      SimplexId downVertex = tree_->getNode(downNodeId)->getVertexId();
      float coordDown[3];
      triangulation_->getVertexPoint(
        downVertex, coordDown[0], coordDown[1], coordDown[2]);

      regionSize = tree_->getNumberOfVisibleRegularNode(i);
      regionSpan = Geometry::distance(coordUp, coordDown);
      regionId = currentZone++;
      // regionId = i;

      // cout << "arc : " << tree_->printArc(i);
      // cout << " span : " << regionSpan;
      // cout << " coords : ";
      // cout << coordDown[0] << ",";
      // cout << coordDown[1] << ",";
      // cout << coordDown[2] << " || ";
      // cout << coordUp[0] << ",";
      // cout << coordUp[1] << ",";
      // cout << coordUp[2] << endl;

      scalarsRegionId->SetTuple1(
        tree_->getNode(downNodeId)->getVertexId(), regionId);
      scalarsRegionId->SetTuple1(
        tree_->getNode(upNodeId)->getVertexId(), regionId);

      scalarsRegionSize->SetTuple1(
        tree_->getNode(downNodeId)->getVertexId(), regionSize);
      scalarsRegionSize->SetTuple1(
        tree_->getNode(upNodeId)->getVertexId(), regionSize);

      scalarsRegionSpan->SetTuple1(
        tree_->getNode(downNodeId)->getVertexId(), regionSpan);
      scalarsRegionSpan->SetTuple1(
        tree_->getNode(upNodeId)->getVertexId(), regionSpan);

      for(SimplexId j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes();
          ++j) {
        SimplexId nodeId = tree_->getSuperArc(i)->getRegularNodeId(j);
        SimplexId vertexId = nodeId;
        // cout << vertexId << ", ";
        if(tree_->getSuperArc(i)->isMasqued(j)) {
          // cout << vertexId << ", ";
          continue;
        }

        // cout << vertexId << ", ";
        scalarsRegionId->SetTuple1(vertexId, regionId);
        scalarsRegionSize->SetTuple1(vertexId, regionSize);
        scalarsRegionSpan->SetTuple1(vertexId, regionSpan);
      }
      // cout << endl;

      // RegionType
      if((upNodeType == CriticalType::Local_minimum
          && downNodeType == CriticalType::Local_maximum)
         || (upNodeType == CriticalType::Local_minimum
             || downNodeType == CriticalType::Local_minimum))
        regionType = static_cast<int>(ArcType::Min_arc);
      else if(upNodeType == CriticalType::Local_maximum
              || downNodeType == CriticalType::Local_maximum)
        regionType = static_cast<int>(ArcType::Max_arc);
      else if(upNodeType == CriticalType::Saddle1
              && downNodeType == CriticalType::Saddle1)
        regionType = static_cast<int>(ArcType::Saddle1_arc);
      else if(upNodeType == CriticalType::Saddle2
              && downNodeType == CriticalType::Saddle2)
        regionType = static_cast<int>(ArcType::Saddle2_arc);
      else
        regionType = static_cast<int>(ArcType::Saddle1_saddle2_arc);

      for(SimplexId j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes();
          ++j) {
        SimplexId nodeId = tree_->getSuperArc(i)->getRegularNodeId(j);
        if(tree_->getSuperArc(i)->isMasqued(j)) {
          // Ignore masqued ones
          continue;
        }
        SimplexId vertexId = nodeId;
        scalarsRegionType->SetTuple1(vertexId, regionType);
      }
    }
  }

  // output
  segmentation_->GetPointData()->AddArray(scalarsRegionId);
  segmentation_->GetPointData()->AddArray(scalarsRegionType);
  segmentation_->GetPointData()->AddArray(scalarsRegionSize);
  segmentation_->GetPointData()->AddArray(scalarsRegionSpan);

  this->printMsg("Topological segmentation built", 1.0, t.getElapsedTime(),
                 this->threadNumber_);
  this->printMsg(std::vector<std::vector<std::string>>{
    {"Region type", std::to_string(scalarsRegionType->GetNumberOfTuples())},
    {"Segmentation Id", std::to_string(scalarsRegionId->GetNumberOfTuples())}});

  toComputeSegmentation_ = false;
}

void ttkContourForests::getTree() {
  // sequential params
  this->preconditionTriangulation(triangulation_);
  this->setVertexScalars(ttkUtils::GetVoidPointer(vtkInputScalars_));
  this->setVertexSoSoffsets(vertexSoSoffsets_);

  this->setTreeType(treeType_);
  // parallel params
  this->setLessPartition(lessPartition_);
  this->setThreadNumber(threadNumber_);
  this->setPartitionNum(partitionNum_);
  // simplification params
  this->setSimplificationMethod(simplificationType_);
  this->setSimplificationThreshold(simplificationThreshold_);
  // build
  ttkVtkTemplateMacro(vtkInputScalars_->GetDataType(),
                      triangulation_->getType(),
                      (this->build<VTK_TT, TTK_TT *>(
                        static_cast<TTK_TT *>(triangulation_->getData()))));

  // what is done is no longer to be done
  toComputeContourTree_ = false;
}

void ttkContourForests::updateTree() {
  // polymorphic tree
  switch(treeType_) {
    case TreeType::Join:
      tree_ = this->getJoinTree();
      break;
    case TreeType::Split:
      tree_ = this->getSplitTree();
      break;
    case TreeType::JoinAndSplit:
      tree_ = this->getJoinTree();
      tree_ = this->getSplitTree();
      break;
    case TreeType::Contour:
      tree_ = this;
      break;
  }

  getCriticalPoints();

  toUpdateTree_ = false;
}

int ttkContourForests::RequestData(vtkInformation *ttkNotUsed(request),
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  vtkWarningMacro(
    "DEPRECATED This plugin will be removed in a future release, please use "
    "FTM instead for contour trees and FTR for Reeb graphs.");

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputSkeletonNodes = vtkPolyData::GetData(outputVector, 0);
  auto outputSkeletonArcs = vtkPolyData::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }

  if(!outputSkeletonNodes || !outputSkeletonArcs || !outputSegmentation) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
#endif

  triangulation_ = ttkAlgorithm::GetTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    this->printErr("Input triangulation is NULL.");
    return -1;
  }
#endif

  varyingMesh_ = false;
  if(triangulation_->isEmpty())
    varyingMesh_ = true;

  // init
  if(varyingMesh_ || !numberOfVertices_) {
    numberOfVertices_ = input->GetNumberOfPoints();
  }

  if(varyingMesh_) {
    segmentation_ = input->NewInstance();
    if(segmentation_ != nullptr) {
      segmentation_->ShallowCopy(input);
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input->GetPointData()) {
    this->printErr("Input has no point data.");
    return -2;
  }
#endif

  // scalars
  vtkInputScalars_ = this->GetInputArrayToProcess(0, input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!vtkInputScalars_) {
    this->printErr("Input scalar is NULL.");
    return -2;
  }
#endif

  varyingDataValues_ = (vtkInputScalars_->GetMTime() > GetMTime());
  if(input->GetPointData()) {

    vertexScalars_.resize(numberOfVertices_);
    for(SimplexId j = 0; j < numberOfVertices_; ++j) {
      vertexScalars_[j] = vtkInputScalars_->GetTuple1(j);
    }
  }

  auto result
    = std::minmax_element(vertexScalars_.begin(), vertexScalars_.end());
  double scalarMin = *result.first;
  double scalarMax = *result.second;
  deltaScalar_ = (scalarMax - scalarMin);

  // offsets
  if(varyingMesh_ || varyingDataValues_ || vertexSoSoffsets_ == nullptr) {

    const auto offsets
      = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);

    if(offsets != nullptr) {
      vertexSoSoffsets_
        = static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsets));
    }

    toUpdateVertexSoSoffsets_ = false;
  }

  if(varyingMesh_ || varyingDataValues_ || !isLoaded_) {
    this->printMsg("Convenient data storage loaded", debug::Priority::DETAIL);
    this->printMsg(
      std::vector<std::vector<std::string>>{
        {"#Tuples", std::to_string(vertexScalars_.size())},
        {"#Vertices", std::to_string(numberOfVertices_)},
        {"Min", std::to_string(scalarMin)},
        {"Max", std::to_string(scalarMax)},
      },
      debug::Priority::DETAIL);
  }

  this->printMsg("Launching computation for field `"
                 + std::string{vtkInputScalars_->GetName()} + "'...");

  isLoaded_ = true;

  if(simplificationType_ == 0) {
    simplificationThreshold_ = simplificationThresholdBuffer_ * deltaScalar_;
  } else if(simplificationType_ == 1) {
    double coord0[3], coord1[3], spanTotal;
    double *bounds = input->GetBounds();
    coord0[0] = bounds[0];
    coord1[0] = bounds[1];
    coord0[1] = bounds[2];
    coord1[1] = bounds[3];
    coord0[2] = bounds[4];
    coord1[2] = bounds[5];
    spanTotal = Geometry::distance(coord0, coord1);
    simplificationThreshold_ = simplificationThresholdBuffer_ * spanTotal;
  } else if(simplificationType_ == 2) {
    simplificationThreshold_
      = simplificationThresholdBuffer_ * triangulation_->getNumberOfVertices();
  }

  // ContourForestsTree //
  if(varyingMesh_ || varyingDataValues_ || toComputeContourTree_) {
    clearTree();
    getTree();
    updateTree();
  }

  // Skeleton //
  if(varyingMesh_ || varyingDataValues_ || toComputeSkeleton_) {
#ifndef TTK_ENABLE_KAMIKAZE
    if(tree_ == nullptr) {
      this->printErr("MergeTree pointer is NULL.");
      return -2;
    }
#endif // TTK_ENABLE_KAMIKAZE
    clearSkeleton();
    getSkeleton();
  }

  // Segmentation //
  if(varyingMesh_ || varyingDataValues_ || toComputeSegmentation_)
    getSegmentation(input);

  // Output //

  // skeleton
  outputSkeletonNodes->ShallowCopy(skeletonNodes_);
  outputSkeletonArcs->ShallowCopy(skeletonArcs_);
  // segmentation
  outputSegmentation->ShallowCopy(segmentation_);

  return 1;
}
