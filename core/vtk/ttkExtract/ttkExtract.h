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
#include <ttkMacros.h>

class TTKEXTRACT_EXPORT ttkExtract : public ttkAlgorithm {

public:
  enum class EXTRACTION_MODE {
    AUTO = -1,
    BLOCKS = 0,
    ROWS = 1,
    GEOMETRY = 2,
    ARRAY_VALUES = 3,
    ARRAYS = 4,
    BLOCK_TUPLES = 5
  };

  enum class VALIDATION_MODE {
    LESS_THEN = 0,
    LESS_EQUAL_THEN = 1,
    EQUAL = 2,
    UNEQUAL = 3,
    GREATER_EQUAL_THEN = 4,
    GREATER_THEN = 5
  };

  enum class CELL_MODE { ALL = 0, ANY = 1, SUB = 2 };

private:
  EXTRACTION_MODE ExtractionMode{EXTRACTION_MODE::AUTO};
  VALIDATION_MODE ValidationMode{VALIDATION_MODE::EQUAL};
  CELL_MODE CellMode{CELL_MODE::ANY};

  int OutputType{-1};
  bool ExtractUniqueValues{true};
  std::string ExpressionString{""};
  int ArrayAttributeType{0};
  std::string OutputArrayName{"Data"};
  int ImageExtent[6]{0, 0, 0, 0, 0, 0};

public:
  ttkSetEnumMacro(ExtractionMode, EXTRACTION_MODE);
  vtkGetEnumMacro(ExtractionMode, EXTRACTION_MODE);

  ttkSetEnumMacro(ValidationMode, VALIDATION_MODE);
  vtkGetEnumMacro(ValidationMode, VALIDATION_MODE);

  ttkSetEnumMacro(CellMode, CELL_MODE);
  vtkGetEnumMacro(CellMode, CELL_MODE);

  vtkSetMacro(OutputType, int);
  vtkGetMacro(OutputType, int);

  vtkSetMacro(ExpressionString, std::string);
  vtkGetMacro(ExpressionString, std::string);

  vtkSetMacro(ExtractUniqueValues, bool);
  vtkGetMacro(ExtractUniqueValues, bool);

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