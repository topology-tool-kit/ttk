/// \ingroup vtk
/// \class ttkForEach
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.01.2019
///
/// \brief TTK VTK-filter that iterates over points, cells, blocks, rows, or
/// field data arrays.
///
/// This filter works in conjunction with the ttkEndFor filter to iterate over
/// points, cells, blocks, rows, or field data arrays
///
/// \param Input vktObject either a vtkMultiBlockDataSet, vtkTable, or
/// vtkDataSet \param Output vktObject one element of the input

#pragma once

// VTK Module
#include <ttkForEachModule.h>

// TTK includes
#include <ttkAlgorithm.h>

class TTKFOREACH_EXPORT ttkForEach : public ttkAlgorithm {

private:
  int Mode{0};
  std::string FieldDataName{""};

public:
  vtkSetMacro(Mode, int);
  vtkGetMacro(Mode, int);

  vtkSetMacro(FieldDataName, std::string);
  vtkGetMacro(FieldDataName, std::string);

  static ttkForEach *New();
  vtkTypeMacro(ttkForEach, ttkAlgorithm)

    protected : ttkForEach();
  ~ttkForEach();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};