/// \ingroup vtk
/// \class ttkLDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief Computes a distance matrix using LDistance between several
/// input datasets with the same number of points
///
/// \sa LDistanceMatrix

#pragma once

// VTK Module
#include <ttkLDistanceMatrixModule.h>

// TTK code includes
#include <LDistanceMatrix.h>
#include <ttkAlgorithm.h>

class TTKLDISTANCEMATRIX_EXPORT ttkLDistanceMatrix
  : public ttkAlgorithm,
    protected ttk::LDistanceMatrix {

public:
  static ttkLDistanceMatrix *New();

  vtkTypeMacro(ttkLDistanceMatrix, ttkAlgorithm);

  vtkSetMacro(DistanceType, const std::string &);
  vtkGetMacro(DistanceType, std::string);

protected:
  ttkLDistanceMatrix();
  ~ttkLDistanceMatrix() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  template <typename T>
  int dispatch(std::vector<std::vector<double>> &distanceMatrix,
               const std::vector<vtkDataSet *> &inputData,
               const size_t nPoints);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
