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
/// The input data arrays needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 for the lowerBoundField, 1 for the upperBoundField
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
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
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/uncertainStartingVortex/">
///   Uncertain Starting Vortex example</a> \n

#pragma once

// VTK Module
#include <ttkMandatoryCriticalPointsModule.h>

// ttk code includes
#include <MandatoryCriticalPoints.h>
#include <ttkAlgorithm.h>

class TTKMANDATORYCRITICALPOINTS_EXPORT ttkMandatoryCriticalPoints
  : public ttkAlgorithm,
    protected ttk::MandatoryCriticalPoints {

public:
  static ttkMandatoryCriticalPoints *New();

  vtkTypeMacro(ttkMandatoryCriticalPoints, ttkAlgorithm);

  void SetSimplificationThreshold(double threshold) {
    if(threshold != simplificationThreshold_) {
      simplificationThreshold_ = threshold;
      simplify_ = true;
      Modified();
    }
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

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  double simplificationThreshold_{0.0};
  bool simplify_{true};

  int outputMinimumComponentId_{0};
  int outputJoinSaddleComponentId_{0};
  int outputSplitSaddleComponentId_{0};
  int outputMaximumComponentId_{0};

  bool outputAllMinimumComponents_{true};
  bool outputAllJoinSaddleComponents_{true};
  bool outputAllSplitSaddleComponents_{true};
  bool outputAllMaximumComponents_{true};

  bool computeMinimumOutput_{true};
  bool computeJoinSaddleOutput_{true};
  bool computeSplitSaddleOutput_{true};
  bool computeMaximumOutput_{true};
};
