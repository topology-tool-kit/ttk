/// \ingroup vtk
/// \class ttkContourForests
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>.
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016
///
///\brief TTK VTK-filter that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
/// This filter takes a scalar field attached as point data to a geometry
/// (either 2D or 3D, either regular grids or triangulations) and computes its
/// contour tree. Several outputs are produced to encode the nodes of the tree
/// (as points in 3D space), the arcs of the tree and the data segmentation.
///
/// \param Input Input scalar field, either 2D or 3D, either regular
/// grid or triangulation (vtkDataSet)
/// \param Output0 Output nodes (vtkPolyData)
/// \param Output1 Output arcs (vtkPolyData)
/// \param Output2 Output data segmentation (vtkUnstructuredGrid)
/// \param Output3 Output persistence diagram (vtkUnstructuredGrid)
/// \param Output4 Output persistence curve (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa ttk::cf::ContourForests

#pragma once

// VTK includes
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkContourForestsModule.h>

// vtk wrapper includes
#include <ttkAlgorithm.h>

// base code includes
#include <Geometry.h>

#include "ContourForests.h"
#include "ContourForestsTree.h"
#include "DeprecatedDataTypes.h"

class TTKCONTOURFORESTS_EXPORT ttkContourForests
  : public ttkAlgorithm,
    protected ttk::cf::ContourForests {
public:
  static ttkContourForests *New();
  vtkTypeMacro(ttkContourForests, ttkAlgorithm);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  void SetForceInputOffsetScalarField(bool onOff);

  void SetTreeType(int tree);

  void ShowMin(bool state);
  void ShowMax(bool state);
  void ShowSaddle1(bool state);
  void ShowSaddle2(bool state);

  void ShowArc(bool state);
  void SetArcResolution(int arcResolution);
  void SetPartitionNumber(int partitionNum);
  void SetLessPartition(bool l);

  void SetSkeletonSmoothing(double skeletonSmooth);

  void SetSimplificationType(int type);

  void SetSimplificationThreshold(double simplificationThreshold);

protected:
  ttkContourForests();

  // VTK Interface //
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  void Modified() override;

  // Base //
  bool isCoincident(float p1[], double p2[]);
  bool isCoincident(double p1[], double p2[]);

  // ContourForestsTree //
  void getTree();
  void updateTree();
  ttk::CriticalType getNodeType(ttk::SimplexId id);
  ttk::CriticalType getNodeType(ttk::SimplexId id,
                                ttk::cf::TreeType type,
                                ttk::cf::MergeTree *tree);
  void getCriticalPoints();
  void clearTree();

  // Skeleton //
  void getSkeleton();
  void clearSkeleton();
  void getSkeletonNodes();
  void getSkeletonArcs();
  int getSkeletonScalars(
    const std::vector<double> &scalars,
    std::vector<std::vector<double>> &skeletonScalars) const;

  // Segmentation //
  void getSegmentation(vtkDataSet *input);
  void clearSegmentation();

  int sample(unsigned int samplingLevel);

  int computeBarycenters();
  void computeSkeleton(unsigned int arcRes);
  void smoothSkeleton(unsigned int skeletonSmoothing);
  void smooth(const ttk::SimplexId idArc, bool order);

private:
  // Base //
  bool isLoaded_{};
  bool lessPartition_{true};
  ttk::cf::MergeTree *tree_{};
  vtkSmartPointer<vtkPolyData> skeletonNodes_{vtkPolyData::New()};
  vtkSmartPointer<vtkPolyData> skeletonArcs_{vtkPolyData::New()};
  vtkSmartPointer<vtkDataSet> segmentation_{};

  // Void //
  vtkNew<vtkUnstructuredGrid> voidUnstructuredGrid_{};
  vtkNew<vtkPolyData> voidPolyData_{};

  // Configuration //
  bool ForceInputOffsetScalarField{false};
  bool varyingMesh_{};
  bool varyingDataValues_{};
  ttk::cf::TreeType treeType_{ttk::cf::TreeType::Contour};
  bool showMin_{true};
  bool showMax_{true};
  bool showSaddle1_{true};
  bool showSaddle2_{true};
  bool showArc_{true};
  unsigned int arcResolution_{20};
  int partitionNum_{-1};
  unsigned int skeletonSmoothing_{15};
  int simplificationType_{};
  double simplificationThreshold_{};
  double simplificationThresholdBuffer_{};

  // Computation handles
  bool toUpdateVertexSoSoffsets_{true};
  bool toComputeContourTree_{true};
  bool toUpdateTree_{true};
  bool toComputeSkeleton_{true};
  bool toComputeSegmentation_{true};

  // Convenient storage //
  vtkDataArray *vtkInputScalars_{};
  double deltaScalar_{};
  ttk::SimplexId numberOfVertices_{};
  ttk::Triangulation *triangulation_{};
  ttk::SimplexId *vertexSoSoffsets_{};
  std::vector<ttk::SimplexId> criticalPoints_{};
  std::vector<double> vertexScalars_{};

  // treeType, SuperArc, several vertices list.
  std::vector<std::vector<std::vector<std::vector<ttk::SimplexId>>>> samples_{};
  std::vector<std::vector<std::vector<std::vector<double>>>> barycenters_{};
};
