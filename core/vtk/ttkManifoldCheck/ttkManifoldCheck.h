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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/manifoldCheck/">Manifold
///   check example</a> \n
#pragma once

// VTK includes -- to adapt

// VTK Module
#include <ttkManifoldCheckModule.h>

// ttk code includes
#include <ManifoldCheck.h>
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKMANIFOLDCHECK_EXPORT ttkManifoldCheck : public ttkAlgorithm,
                                                 protected ttk::ManifoldCheck {

public:
  static ttkManifoldCheck *New();
  vtkTypeMacro(ttkManifoldCheck, ttkAlgorithm);

protected:
  ttkManifoldCheck();

  ~ttkManifoldCheck() override{};

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<ttk::SimplexId> vertexLinkComponentNumber_;
  std::vector<ttk::SimplexId> edgeLinkComponentNumber_;
  std::vector<ttk::SimplexId> triangleLinkComponentNumber_;
};
