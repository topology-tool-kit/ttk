/// \ingroup vtk
/// \class ttkManifoldCheck
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2018.
///
/// \brief TTK VTK-filter for manifold checks.
///
/// This filter performs a manifold check for each simplex, by counting the
/// number of connected components of link. On a d-dimensional triangulation,
/// this number should be equal to 1 for all but (d-1)-simplices, for which it
/// can be 1 (boundary simplices) or 2 (interior simplices). Other values
/// indicate a non-manifold simplex.
///
/// The link component number is stored as an integer array for each type of
/// simplex. In practice, these arrays are both stored as (i) point data and
/// (ii) cell data, by considering the maximum value achieved by a co-face (i)
/// and by a vertex (ii) respectively.
///
/// \param Input Input data set (vtkDataSet)
/// \param Output Output  field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::ManifoldCheck
#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkGenericCell.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <ManifoldCheck.h>
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkManifoldCheck
#else
class ttkManifoldCheck
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkManifoldCheck *New();
  vtkTypeMacro(ttkManifoldCheck, vtkDataSetAlgorithm)

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

protected:
  ttkManifoldCheck() {

    // init
    UseAllCores = true;
  }

  ~ttkManifoldCheck(){};

  TTK_SETUP();

private:
  std::vector<ttk::SimplexId> vertexLinkComponentNumber_;
  std::vector<ttk::SimplexId> edgeLinkComponentNumber_;
  std::vector<ttk::SimplexId> triangleLinkComponentNumber_;
  ttk::ManifoldCheck manifoldCheck_;
};
