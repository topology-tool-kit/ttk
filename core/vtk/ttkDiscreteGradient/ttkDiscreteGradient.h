/// \ingroup vtk
/// \class ttkDiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@sorbonne-universite.fr>
/// \date April 2018.
///
/// \brief TTK VTK-filter that wraps the discreteGradient processing package.
///
/// VTK wrapping code for the @DiscreteGradient package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::DiscreteGradient
#ifndef _TTK_DISCRETEGRADIENT_H
#define _TTK_DISCRETEGRADIENT_H

// VTK includes
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
#include <vtkLine.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

// VTK Module
#include <ttkDiscreteGradientModule.h>

// ttk code includes
#include <DiscreteGradient.h>
#include <ttkTriangulationAlgorithm.h>

class TTKDISCRETEGRADIENT_EXPORT ttkDiscreteGradient
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkDiscreteGradient *New();

  vtkTypeMacro(ttkDiscreteGradient, vtkDataSetAlgorithm);

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

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(ComputeGradientGlyphs, bool);
  vtkGetMacro(ComputeGradientGlyphs, bool);

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  int setupTriangulation(vtkDataSet *input);
  int getScalars(vtkDataSet *input);
  int getOffsets(vtkDataSet *input);

protected:
  ttkDiscreteGradient();
  ~ttkDiscreteGradient() override;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  template <typename T>
  int dispatch(vtkUnstructuredGrid *outputCriticalPoints);

  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  bool ForceInputOffsetScalarField;
  bool ComputeGradientGlyphs;
  int IterationThreshold;
  int ScalarFieldId;
  int OffsetFieldId;

  ttk::Triangulation *triangulation_;
  ttk::dcg::DiscreteGradient discreteGradient_;
  vtkDataArray *inputScalars_;
  ttkSimplexIdTypeArray *offsets_;
  vtkDataArray *inputOffsets_;
  bool hasUpdatedMesh_;
};

#endif // _TTK_DISCRETEGRADIENT_H
