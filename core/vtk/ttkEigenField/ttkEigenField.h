/// \ingroup vtk
/// \class ttkEigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2019
///
/// \brief TTK VTK-filter for harmonic field computations.
///
/// The current filter takes a list of sources with attached scalar
/// values and produces a scalar harmonic field fulfilling these
/// constraints.
///
/// \param Input0 Input geometry, either 2D or 3D, either regular grid
/// or triangulation (vtkDataSet)
/// \param Input1 List of critical point constraints (vtkPointSet)
/// \param Output Output harmonic scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Dynamic harmonic fields for surface processing" \n
/// Kai Xu, Hao Zhang, Daniel Cohen-Or, Yueshan Xiong \n
/// Computers & Graphics 2009. \n
///
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkFTMTree
/// \sa ttkIdentifiers
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

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(EigenNumber, size_t);
  vtkGetMacro(EigenNumber, size_t);

  // get mesh from VTK
  int getTriangulation(vtkDataSet *input);
  // get array of identifiers on the mesh
  int getIdentifiers(vtkPointSet *input);
  // get constraint values on identifiers
  int getConstraints(vtkPointSet *input);

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
  // output scalar field
  std::string OutputScalarFieldName{"OutputEigenField"};
  // number of eigenpairs to compute
  size_t EigenNumber{20};

  // enum: float or double
  int OutputScalarFieldType{EigenFieldType::Float};
  // worker object
  ttk::EigenField baseWorker_{};
  // teh mesh
  ttk::Triangulation *triangulation_{};
};
