/// \ingroup vtk
/// \class ttkEigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date April 2019
///
/// \brief TTK VTK-filter for eigenfunctions computation.
///
/// This plugin computes the first eigenfunctions of a given
/// triangular surface mesh.
///
/// \param Input0 Input 2D geometry, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output eigenfunctions (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
///  "Spectral surface quadrangulation"
///  Shen Dong, Peer-Timo Bremer, Michael Garland, Valerio Pascucci, John C.
///  Hart SIGGRAPH 2006
///
/// \sa ttkHarmonicField
/// \sa ttk::EigenField

#pragma once

// VTK includes -- to adapt
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

// ttk code includes
#include <EigenField.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

enum EigenFieldType { Float = 0, Double };

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkEigenField
#else
class ttkEigenField
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkEigenField *New();
  vtkTypeMacro(ttkEigenField, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(OutputFieldName, std::string);
  vtkGetMacro(OutputFieldName, std::string);

  vtkSetMacro(EigenNumber, unsigned int);
  vtkGetMacro(EigenNumber, unsigned int);

  vtkSetMacro(ComputeStatistics, bool);
  vtkGetMacro(ComputeStatistics, bool);

  // get mesh from VTK
  int getTriangulation(vtkDataSet *input);

  // default copy constructor
  ttkEigenField(const ttkEigenField &) = delete;
  // default move constructor
  ttkEigenField(ttkEigenField &&) = delete;
  // default copy assignment operator
  ttkEigenField &operator=(const ttkEigenField &) = delete;
  // default move assignment operator
  ttkEigenField &operator=(ttkEigenField &&) = delete;

protected:
  ttkEigenField();

  ~ttkEigenField() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  // output field name
  std::string OutputFieldName{"OutputEigenFunctions"};
  // number of eigenpairs to compute
  unsigned int EigenNumber{500};
  // if statistics are to be computed
  bool ComputeStatistics{false};

  // enum: float or double
  int OutputFieldType{EigenFieldType::Float};
  // worker object
  ttk::EigenField baseWorker_{};
  // teh mesh
  ttk::Triangulation *triangulation_{};
};
