/// \ingroup vtk
/// \class ttkCompactTriangulationPreconditioning
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date May 2021.
///
/// \brief TTK VTK-filter for the preconditioning of the compact triangulation.
///
/// Given a simplicial mesh, this filter uses the PR star octree to divide
/// the mesh into different regions, and adds this clustering information as
/// a new scalar field to the original dataset. This clustering index scalar
/// field can be further used by Compact Triangulation.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// See the corresponding standalone program for a usage example:
///   - standalone/CompactTriangulationPreconditioning/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::CompactTriangulationPreconditioning
/// \sa ttk::TopoCluster

#pragma once

// VTK Module
#include <ttkCompactTriangulationPreconditioningModule.h>

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

// TTK Base Includes
#include <CompactTriangulationPreconditioning.h>

class TTKCOMPACTTRIANGULATIONPRECONDITIONING_EXPORT
  ttkCompactTriangulationPreconditioning
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::CompactTriangulationPreconditioning // and we inherit from
                                                       // the base class
{
private:
  int Threshold;
  std::vector<std::string> scalarFields;

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
  static ttkCompactTriangulationPreconditioning *New();
  vtkTypeMacro(ttkCompactTriangulationPreconditioning, ttkAlgorithm);

protected:
  ttkCompactTriangulationPreconditioning();
  ~ttkCompactTriangulationPreconditioning() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};