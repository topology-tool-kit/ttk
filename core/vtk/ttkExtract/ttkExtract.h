/// \ingroup vtk
/// \class ttkExtract
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that extracts blocks or a geometry subset
///
/// This filter uses a list of values to extract either blocks of a
/// 'vtkMultiBlockDataSet' by interpreting the values as block indices, or the
/// subset of a 'vtkDataObject' whose point/cell values are contained in that
/// list.
///
#pragma once

// VTK Module
#include <ttkExtractModule.h>

// TTK includes
#include <ttkAlgorithm.h>

class TTKEXTRACT_EXPORT ttkExtract : public ttkAlgorithm {

private:
  int Mode{0};
  int OutputType{0};
  std::string ExpressionString{""};
  int CellMode{0};
  double ImageBounds[6]{0, 0, 0, 0, 0, 0};

public:
  vtkSetMacro(Mode, int);
  vtkGetMacro(Mode, int);

  vtkSetMacro(OutputType, int);
  vtkGetMacro(OutputType, int);

  vtkSetMacro(ExpressionString, std::string);
  vtkGetMacro(ExpressionString, std::string);

  vtkSetMacro(CellMode, int);
  vtkGetMacro(CellMode, int);

  vtkSetVector6Macro(ImageBounds, double);
  vtkGetVector6Macro(ImageBounds, double);

  int GetVtkDataTypeName(int outputType, std::string &dataTypeName);

  static ttkExtract *New();
  vtkTypeMacro(ttkExtract, ttkAlgorithm)

    protected : ttkExtract();
  ~ttkExtract();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};