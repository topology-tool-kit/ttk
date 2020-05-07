/// \ingroup vtk
/// \class ttkBarycentricSubdivision
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the barycentricSubdivision processing
/// package.
///
/// VTK wrapping code for the @BarycentricSubdivision package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::BarycentricSubdivision
#pragma once

// VTK Module
#include <ttkBarycentricSubdivisionModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <BarycentricSubdivision.h>

class TTKBARYCENTRICSUBDIVISION_EXPORT ttkBarycentricSubdivision
  : public ttkAlgorithm,
    protected ttk::BarycentricSubdivision {

private:
  // number of subdivisions
  unsigned int SubdivisionLevel{1};

  // output 3D coordinates of generated points: old points first, then edge
  // middles, then triangle barycenters
  std::vector<float> points_{};

  // output triangles
  std::vector<ttk::LongSimplexId> cells_{};

  // generated point cell id
  std::vector<ttk::SimplexId> pointId_{};

  // generated points dimension: 0 vertex of parent triangulation, 1 edge
  // middle, 2 triangle barycenter
  std::vector<ttk::SimplexId> pointDim_{};

public:
  vtkGetMacro(SubdivisionLevel, unsigned int);
  vtkSetMacro(SubdivisionLevel, unsigned int);

  static ttkBarycentricSubdivision *New();
  vtkTypeMacro(ttkBarycentricSubdivision, ttkAlgorithm);

  vtkSmartPointer<vtkDataArray>
    AllocateScalarField(vtkDataArray *const inputScalarField,
                        int ntuples) const;

  int InterpolateScalarFields(vtkUnstructuredGrid *const input,
                              vtkUnstructuredGrid *const output) const;

protected:
  ttkBarycentricSubdivision();
  ~ttkBarycentricSubdivision();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
