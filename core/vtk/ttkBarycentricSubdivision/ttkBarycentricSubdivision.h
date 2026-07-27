/// \ingroup vtk
/// \class ttkBarycentricSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2019.
///
/// \brief TTK VTK-filter that wraps the ttk::BarycentricSubdivision
/// processing package.
///
/// VTK wrapping code for the ttk::BarycentricSubdivision package.
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

// TTK code includes
#include <BarycentricSubdivision.h>
#include <ttkAlgorithm.h>

class vtkUnstructuredGrid;

template <typename T>
class vtkSmartPointer;

class TTKBARYCENTRICSUBDIVISION_EXPORT ttkBarycentricSubdivision
  : public ttkAlgorithm,
    protected ttk::BarycentricSubdivision {

public:
  static ttkBarycentricSubdivision *New();
  vtkGetMacro(SubdivisionLevel, unsigned int);
  vtkSetMacro(SubdivisionLevel, unsigned int);

  vtkTypeMacro(ttkBarycentricSubdivision, ttkAlgorithm);

protected:
  ttkBarycentricSubdivision();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  /**
   * @brief Allocate an output array of same type that input array
   */
  vtkSmartPointer<vtkDataArray>
    AllocateScalarField(vtkDataArray *const inputScalarField,
                        int ntuples) const;

  int InterpolateScalarFields(vtkDataSet *const input,
                              vtkUnstructuredGrid *const output,
                              ttk::Triangulation &inputTriangulation) const;

private:
  // number of subdivisions
  unsigned int SubdivisionLevel{1};
};
