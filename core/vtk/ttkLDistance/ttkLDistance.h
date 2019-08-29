/// \ingroup vtk
/// \class ttkLDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK VTK-filter that wraps the lDistance processing package.
///
/// VTK wrapping code for the @LDistance package.
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
/// \sa ttk::LDistance
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
#include <vtkSmartPointer.h>

// ttk code includes
#include <LDistance.h>
#include <ttkWrapper.h>

// In this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkLDistance
#else
class ttkLDistance
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkLDistance *New();

  vtkTypeMacro(ttkLDistance, vtkDataSetAlgorithm);

  // Default ttk setters
  vtkSetMacro(debugLevel_, int);

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
  vtkSetMacro(ScalarField1, std::string);
  vtkGetMacro(ScalarField1, std::string);

  vtkSetMacro(ScalarField2, std::string);
  vtkGetMacro(ScalarField2, std::string);

  vtkSetMacro(ScalarFieldId1, int);
  vtkGetMacro(ScalarFieldId1, int);

  vtkSetMacro(ScalarFieldId2, int);
  vtkGetMacro(ScalarFieldId2, int);

  vtkSetMacro(DistanceType, std::string);
  vtkGetMacro(DistanceType, std::string);

  vtkSetMacro(DistanceFieldName, std::string);
  vtkGetMacro(DistanceFieldName, std::string);

  vtkGetMacro(result, double);

protected:
  // By default, this filter has one input and one output.
  ttkLDistance() {

    DistanceType = "2";
    DistanceFieldName = "L2-distance";
    ScalarFieldId1 = 0;
    ScalarFieldId2 = 1;
    outputScalarField_ = NULL;
    UseAllCores = true;
    result = -1.;
  }

  ~ttkLDistance(){};

  TTK_SETUP();

private:
  std::string DistanceType;
  std::string ScalarField1;
  std::string ScalarField2;
  int ScalarFieldId1;
  int ScalarFieldId2;
  std::string DistanceFieldName;
  double result;

  vtkSmartPointer<vtkDataArray> outputScalarField_;
  ttk::LDistance lDistance_;
};
