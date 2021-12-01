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
/// \b Related \b publications \n
/// "The PR-star octree: A spatio-topological data structure for tetrahedral
/// meshes." Kenneth Weiss, Leila Floriani, Riccardo Fellegara, and Marcelo
/// Velloso Proc. of ACM SIGSPATIAL 2011.
///
/// "TopoCluster: A Localized Data Structure for Topology-based Visualization"
/// Guoxi Liu, Federico Iuricich, Riccardo Fellegara, and Leila De Floriani
/// IEEE Transactions on Visualization and Computer Graphics, 2021.
///
/// \sa ttk::CompactTriangulationPreconditioning
/// \sa ttk::TopoCluster

#pragma once

// VTK Module
#include <ttkCompactTriangulationPreconditioningModule.h>

// VTK Includes
// VTK includes -- to adapt
#include <ttkAlgorithm.h>
#include <vtkDataArraySelection.h>
#include <vtkSmartPointer.h>

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
  vtkSmartPointer<vtkDataArraySelection> ArraySelection;

public:
  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkCompactTriangulationPreconditioning *New();
  vtkTypeMacro(ttkCompactTriangulationPreconditioning, ttkAlgorithm);

  vtkSetMacro(Threshold, int);
  vtkGetMacro(Threshold, int);

  // copy the vtkPassSelectedArray ("PassArrays" filter) API
  vtkDataArraySelection *GetDataArraySelection() {
    return this->ArraySelection.GetPointer();
  }

  void SetDataArraySelection(
    const vtkSmartPointer<vtkDataArraySelection> &selection) {
    this->ArraySelection = selection;
  }

protected:
  ttkCompactTriangulationPreconditioning();
  ~ttkCompactTriangulationPreconditioning() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};