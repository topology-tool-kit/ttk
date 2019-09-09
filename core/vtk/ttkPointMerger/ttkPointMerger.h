/// \ingroup vtk
/// \class ttkPointMerger
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2018.
///
/// \brief TTK VTK-filter for point merging.
///
/// This filter merges the points of a mesh whose distance is lower than a user
/// defined threshold.
///
/// \param Input Input data set (vtkDataSet)
/// \param Output Output data set (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <Geometry.h>
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPointMerger
#else
class ttkPointMerger
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkPointMerger *New();
  vtkTypeMacro(ttkPointMerger, vtkDataSetAlgorithm)

    // default ttk setters
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

  vtkSetMacro(BoundaryOnly, bool);
  vtkGetMacro(BoundaryOnly, bool);

  vtkSetMacro(DistanceThreshold, double);
  vtkGetMacro(DistanceThreshold, double);

protected:
  ttkPointMerger() {

    // init
    DistanceThreshold = 0.001;
    BoundaryOnly = true;
    UseAllCores = true;
  }

  ~ttkPointMerger(){};

  TTK_SETUP();

private:
  bool BoundaryOnly;
  double DistanceThreshold;
};
