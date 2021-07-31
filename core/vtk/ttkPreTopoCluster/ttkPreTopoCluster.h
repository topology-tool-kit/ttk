/// \ingroup vtk
/// \class ttkPreTopoCluster
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date May 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::PreTopoCluster module.
///
/// Given a simplicial mesh, this filter uses the PR star octree to divide
/// the mesh into different regions, and adds this clustering information as
/// a new scalar field to the original dataset. This clustering index scalar
/// field can be further used by TopoCluster data structure.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// See the corresponding standalone program for a usage example:
///   - standalone/PreTopoCluster/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::PreTopoCluster
/// \sa ttk::TopoCluster

#pragma once

// VTK Module
#include <ttkPreTopoClusterModule.h>

// VTK Includes
// VTK includes -- to adapt
#include <ttkAlgorithm.h>
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
#include <vtkUnstructuredGrid.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkPreTopoCluster
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <PreTopoCluster.h>

class TTKPRETOPOCLUSTER_EXPORT ttkPreTopoCluster
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::PreTopoCluster // and we inherit from the base class
{
private:
  int Threshold;
  vector<string> scalarFields;

public:
  vtkSetMacro(Threshold, int);
  vtkGetMacro(Threshold, int);

  void SetScalarField(string name) {
    scalarFields.push_back(name);
    Modified();
  }

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkPreTopoCluster *New();
  vtkTypeMacro(ttkPreTopoCluster, ttkAlgorithm);

protected:
  ttkPreTopoCluster();
  ~ttkPreTopoCluster() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};