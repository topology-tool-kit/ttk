/// \ingroup vtk
/// \class ttkLDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief Computes a distance matrix using LDistance between several
/// same-dimensions input datasets
///
/// \sa LDistanceMatrix

#pragma once

// VTK includes -- to adapt
#include <vtkInformation.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

// VTK Module
#include <ttkLDistanceMatrixModule.h>

// TTK code includes
#include <LDistanceMatrix.h>
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKLDISTANCEMATRIX_EXPORT ttkLDistanceMatrix
  : public vtkMultiBlockDataSetAlgorithm,
    protected ttk::LDistanceMatrix {

public:
  static ttkLDistanceMatrix *New();

  vtkTypeMacro(ttkLDistanceMatrix, vtkMultiBlockDataSetAlgorithm);

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(DistanceType, std::string);
  vtkGetMacro(DistanceType, std::string);

protected:
  ttkLDistanceMatrix();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string ScalarField{};
};
