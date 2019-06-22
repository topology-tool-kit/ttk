/// \ingroup vtk
/// \class ttkMandatoryCriticalPoints
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2016.
///
/// \brief TTK VTK-filter for the computation of mandatory critical
/// points in uncertain scalar data.
///
/// This filter computes the mandatory critical points of uncertain scalar
/// fields defined on triangulations. The input uncertain data is represented
/// by reliable bound fields for each vertex. In particular, the input geometry
/// must be associated with two point data scalar fields, representing the
/// lower and upper bounds for each vertex.
///
/// The output is a domain segmentation into the mandatory critical points.
///
/// \param Input Input uncertain scalar field represented by lower and upper
/// bounds, either 2D or 3D, either regular grid or triangulation (vtkDataSet)
/// triangulation (vtkDataSet)
/// \param Output0 Output mandatory minimum segmentation (vtkDataSet)
/// \param Output1 Output mandatory join saddle segmentation (vtkDataSet)
/// \param Output2 Output mandatory split saddle segmentation (vtkDataSet)
/// \param Output3 Output mandatory maximum segmentation (vtkDataSet)
/// \param Output4 Output mandatory join tree (vtkUnstructuredGrid)
/// \param Output5 Output mandatory split tree (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Mandatory Critical Points of 2D Uncertain Scalar Fields" \n
/// David Guenther, Joseph Salmon, Julien Tierny \n
/// Proc. of EuroVis 2014. \n
/// Computer Graphics Forum, 2014.
///
/// \sa ttk::MandatoryCriticalPoints
/// \sa vtkUncertainDataEstimator
///
#ifndef _TTK_MANDATORYCRITICALPOINTS_H
#define _TTK_MANDATORYCRITICALPOINTS_H

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <MandatoryCriticalPoints.h>
#include <ttkWrapper.h>

#include <queue>
#include <utility>
#include <vector>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMandatoryCriticalPoints
#else
class ttkMandatoryCriticalPoints
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkMandatoryCriticalPoints *New();

  vtkTypeMacro(ttkMandatoryCriticalPoints, vtkDataSetAlgorithm);

  // int buildMandatoryTree(vtkUnstructuredGrid *outputJoinTree, const Digraph
  // &mandatoryTree, const bool isJoinTree); void
  // buildMandatorySplitTree(vtkUnstructuredGrid *outputSplitTree);

  int buildVtkTree(vtkUnstructuredGrid *outputTree,
                   ttk::MandatoryCriticalPoints::TreeType treeType);

  void SetDebugLevel(int debugLevel) {
    if(debugLevel != debugLevel_) {
      debugLevel_ = debugLevel;
      computeAll_ = true;
      Modified();
    }
  }

  void SetSimplificationThreshold(double threshold) {
    if(threshold != simplificationThreshold_) {
      simplificationThreshold_ = threshold;
      simplify_ = true;
      Modified();
    }
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
    computeAll_ = true;
  }

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.

  void SetUpperBoundField(const int &id) {
    upperBoundId = id;
    Modified();
  }

  void SetUpperBoundFieldName(std::string name) {
    upperBoundFiledName_ = name;
    Modified();
  }

  void SetLowerBoundField(const int &id) {
    lowerBoundId = id;
    Modified();
  }

  void SetLowerBoundFieldName(std::string name) {
    lowerBoundFieldName_ = name;
    Modified();
  }

  void SetOutputMinimumComponentId(int id) {
    outputMinimumComponentId_ = id;
    computeMinimumOutput_ = true;
    Modified();
  }

  void SetOutputJoinSaddleComponentId(int id) {
    outputJoinSaddleComponentId_ = id;
    computeJoinSaddleOutput_ = true;
    Modified();
  }

  void SetOutputSplitSaddleComponentId(int id) {
    outputSplitSaddleComponentId_ = id;
    computeSplitSaddleOutput_ = true;
    Modified();
  }

  void SetOutputMaximumComponentId(int id) {
    outputMaximumComponentId_ = id;
    computeMaximumOutput_ = true;
    Modified();
  }

  void setOutputAllMinimumComponents(bool outputAll) {
    outputAllMinimumComponents_ = outputAll;
    computeMinimumOutput_ = true;
    Modified();
  }

  void setOutputAllJoinSaddleComponents(bool outputAll) {
    outputAllJoinSaddleComponents_ = outputAll;
    computeJoinSaddleOutput_ = true;
    Modified();
  }

  void setOutputAllSplitSaddleComponents(bool outputAll) {
    outputAllSplitSaddleComponents_ = outputAll;
    computeSplitSaddleOutput_ = true;
    Modified();
  }

  void setOutputAllMaximumComponents(bool outputAll) {
    outputAllMaximumComponents_ = outputAll;
    computeMaximumOutput_ = true;
    Modified();
  }

protected:
  ttkMandatoryCriticalPoints();

  ~ttkMandatoryCriticalPoints();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  TTK_PIPELINE_REQUEST();
  TTK_OUTPUT_MANAGEMENT();

private:
  ttk::MandatoryCriticalPoints mandatoryCriticalPoints_;

  unsigned long int inputMTime_;
  bool computeAll_;

  bool UseAllCores;
  int ThreadNumber, lowerBoundId, upperBoundId;
  std::string upperBoundFiledName_;
  std::string lowerBoundFieldName_;
  vtkIntArray *outputMandatoryMinimum_;
  vtkIntArray *outputMandatoryJoinSaddle_;
  vtkIntArray *outputMandatorySplitSaddle_;
  vtkIntArray *outputMandatoryMaximum_;

  // Join Trees
  vtkPoints *mandatoryJoinTreePoints_;
  std::vector<vtkIdList *> mandatoryJoinTreeEdge_;
  vtkIntArray *mdtJoinTreePointType_;
  vtkDoubleArray *mdtJoinTreePointLowInterval_;
  vtkDoubleArray *mdtJoinTreePointUpInterval_;
  vtkIntArray *mdtJoinTreePointComponentId_;
  vtkIntArray *mdtJoinTreeEdgeSwitchable_;
  // Split Tree
  vtkPoints *mandatorySplitTreePoints_;
  std::vector<vtkIdList *> mandatorySplitTreeEdge_;
  vtkIntArray *mdtSplitTreePointType_;
  vtkDoubleArray *mdtSplitTreePointLowInterval_;
  vtkDoubleArray *mdtSplitTreePointUpInterval_;
  vtkIntArray *mdtSplitTreePointComponentId_;
  vtkIntArray *mdtSplitTreeEdgeSwitchable_;

  double simplificationThreshold_;
  bool simplify_;

  int outputMinimumComponentId_;
  int outputJoinSaddleComponentId_;
  int outputSplitSaddleComponentId_;
  int outputMaximumComponentId_;

  bool outputAllMinimumComponents_;
  bool outputAllJoinSaddleComponents_;
  bool outputAllSplitSaddleComponents_;
  bool outputAllMaximumComponents_;

  bool computeMinimumOutput_;
  bool computeJoinSaddleOutput_;
  bool computeSplitSaddleOutput_;
  bool computeMaximumOutput_;

  ttk::Triangulation *triangulation_;

  float memoryUsage_;

  // base code features
  int doIt(std::vector<vtkDataSet *> &inputs,
           std::vector<vtkDataSet *> &outputs);

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

#endif // _TTK_MANDATORYCRITICALPOINTS_H
