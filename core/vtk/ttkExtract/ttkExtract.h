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
  int ExtractionMode{0};
  int OutputType{-1};
  bool ExtractUniqueValues{true};
  std::string ExpressionString{""};
  int ValidationMode{0};
  int CellMode{0};
  int ArrayAttributeType{0};
  std::string OutputArrayName{"Data"};
  int ImageExtent[6]{0, 0, 0, 0, 0, 0};

public:
  vtkSetMacro(ExtractionMode, int);
  vtkGetMacro(ExtractionMode, int);

  vtkSetMacro(OutputType, int);
  vtkGetMacro(OutputType, int);

  vtkSetMacro(ExpressionString, std::string);
  vtkGetMacro(ExpressionString, std::string);

  vtkSetMacro(ExtractUniqueValues, bool);
  vtkGetMacro(ExtractUniqueValues, bool);

  vtkSetMacro(ValidationMode, int);
  vtkGetMacro(ValidationMode, int);

  vtkSetMacro(CellMode, int);
  vtkGetMacro(CellMode, int);

  vtkSetMacro(ArrayAttributeType, int);
  vtkGetMacro(ArrayAttributeType, int);

  vtkSetMacro(OutputArrayName, std::string);
  vtkGetMacro(OutputArrayName, std::string);

  vtkSetVector6Macro(ImageExtent, int);
  vtkGetVector6Macro(ImageExtent, int);

  std::string GetVtkDataTypeName(const int outputType) const;

  int ExtractBlocks(vtkDataObject *output,
                    vtkDataObject *input,
                    const std::vector<double> &indices,
                    const bool &extractTuples) const;

  int ExtractRows(vtkDataObject *output,
                  vtkDataObject *input,
                  const std::vector<double> &indices) const;

  int ExtractGeometry(vtkDataObject *output,
                      vtkDataObject *input,
                      const std::vector<double> &labels);

  int ExtractArrayValues(vtkDataObject *output,
                         vtkDataObject *input,
                         const std::vector<double> &indices);

  int ExtractArray(vtkDataObject *output,
                   vtkDataObject *input,
                   const std::vector<double> &indices);

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