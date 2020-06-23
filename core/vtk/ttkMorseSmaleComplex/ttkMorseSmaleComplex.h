/// \ingroup vtk
/// \class ttkMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK VTK-filter that wraps the morseSmaleComplex processing package.
///
/// TTK module for the computation of Morse-Smale complexes.
/// Morse-Smale complexes are useful topological abstractions of scalar
/// fields for data segmentation, feature extraction, etc.
///
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkUnstructuredGrid)
/// \param Output1 Output 1-separatrices (vtkUnstructuredGrid)
/// \param Output2 Output 2-separatrices (vtkUnstructuredGrid)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MorseSmaleComplex
///
#ifndef _TTK_MORSESMALECOMPLEX_H
#define _TTK_MORSESMALECOMPLEX_H

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkMorseSmaleComplexModule.h>

// ttk code includes
#include <MorseSmaleComplex.h>
#include <ttkTriangulationAlgorithm.h>

class TTKMORSESMALECOMPLEX_EXPORT ttkMorseSmaleComplex
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkMorseSmaleComplex *New();

  vtkTypeMacro(ttkMorseSmaleComplex, vtkDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
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

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(PeriodicBoundaryConditions, int);
  vtkGetMacro(PeriodicBoundaryConditions, int);

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(ComputeCriticalPoints, bool);
  vtkGetMacro(ComputeCriticalPoints, bool);

  vtkSetMacro(ComputeAscendingSeparatrices1, bool);
  vtkGetMacro(ComputeAscendingSeparatrices1, bool);

  vtkSetMacro(ComputeDescendingSeparatrices1, bool);
  vtkGetMacro(ComputeDescendingSeparatrices1, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkSetMacro(ComputeAscendingSeparatrices2, bool);
  vtkGetMacro(ComputeAscendingSeparatrices2, bool);

  vtkSetMacro(ComputeDescendingSeparatrices2, bool);
  vtkGetMacro(ComputeDescendingSeparatrices2, bool);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

  vtkSetMacro(ReturnSaddleConnectors, int);
  vtkGetMacro(ReturnSaddleConnectors, int);

  vtkSetMacro(SaddleConnectorsPersistenceThreshold, double);
  vtkGetMacro(SaddleConnectorsPersistenceThreshold, double);

  int setupTriangulation(vtkDataSet *input);
  vtkDataArray *getScalars(vtkDataSet *input);
  vtkDataArray *getOffsets(vtkDataSet *input);

protected:
  template <typename VTK_TT>
  int dispatch(vtkDataArray *inputScalars,
               vtkDataArray *inputOffsets,
               vtkUnstructuredGrid *outputCriticalPoints,
               vtkUnstructuredGrid *outputSeparatrices1,
               vtkUnstructuredGrid *outputSeparatrices2);

  ttkMorseSmaleComplex();
  ~ttkMorseSmaleComplex() override;

  TTK_SETUP();

  virtual int FillInputPortInformation(int port, vtkInformation *info) override;
  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override;

private:
  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  bool ForceInputOffsetScalarField;
  bool PeriodicBoundaryConditions;
  int IterationThreshold;
  bool ComputeCriticalPoints;
  bool ComputeAscendingSeparatrices1;
  bool ComputeDescendingSeparatrices1;
  bool ComputeSaddleConnectors;
  bool ComputeAscendingSeparatrices2;
  bool ComputeDescendingSeparatrices2;
  bool ComputeAscendingSegmentation;
  bool ComputeDescendingSegmentation;
  bool ComputeFinalSegmentation;
  int ScalarFieldId;
  int OffsetFieldId;
  int ReturnSaddleConnectors;
  double SaddleConnectorsPersistenceThreshold;

  ttk::MorseSmaleComplex morseSmaleComplex_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *defaultOffsets_;
  bool hasUpdatedMesh_;
};

#endif // _TTK_MORSESMALECOMPLEX_H
