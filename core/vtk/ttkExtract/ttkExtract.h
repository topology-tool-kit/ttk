/// \ingroup vtk
/// \class ttkExtract
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that provides multiple methods to extract subsets of
/// an input data object based on a logical expression.
///
/// This filter provides multiple methods to extract subsets of an input data
/// object based on a logical expression:
///
/// 1. Blocks: The filter extracts all blocks of a vtkMultiBlockDataSet based on
/// a list of block indices. It is also possible to extract a single block of a
/// vtkMultiBlockDataSet and explicitly specify its type, which is then returned
/// instead of a vtkMultiBlockDataSet containing a single block. This is
/// especially useful to extract vtkImageData objects.
///
/// 2. Block Tuples: Many pipelines produce vtkMultiBlockDataSets that contain
/// vtkMultiBlockDataSets that represent lists. For example, a parent
/// vtkMultiBlockDataSet might contain lists of Merge Trees, Persistence
/// Diagrams, and Domain segmentations, where each entry in a list represents a
/// timestep/ensemble member. Extracting all elements for a given list of
/// timesteps/ensemble members is very cumbersome with the original block
/// extraction method. The block tuples mode makes it possible to conveniently
/// extract these tuples based on a list of timesteps/ensemble member indices.
///
/// 3. Rows: The filter extracts all rows of a vtkTable based on a list of row
/// indices.
///
/// 4. Geometry: The filter extracts the subset of the input geometry whose
/// point/cell data satisfies a logical expression. It is also possible to pass
/// on the input dataset and only add a mask array that marks points/cells that
/// satisfy the condition.
///
/// 5. Array Values: The filter extracts all array values of a vtkAbstractArray
/// based on a list of value indices. The extracted values are stored in a new
/// field data array.
///
/// 6. Arrays: The filter extracts all point/cell data arrays based on a given
/// list of indices (not names).
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

  enum class CELL_MODE { ALL = 0, ANY = 1 };

private:
  EXTRACTION_MODE ExtractionMode{EXTRACTION_MODE::AUTO};
  VALIDATION_MODE ValidationMode{VALIDATION_MODE::EQUAL};
  CELL_MODE CellMode{CELL_MODE::ALL};

  int OutputType{-1};
  std::string ExpressionString{""};
  int ArrayAttributeType{0};
  bool ExtractUniqueValues{true};
  std::string OutputArrayName{"Data"};
  int ImageExtent[6]{0, 0, 0, 0, 0, 0};
  bool MaskOnly{false};

public:
  ttkSetEnumMacro(ExtractionMode, EXTRACTION_MODE);
  vtkGetEnumMacro(ExtractionMode, EXTRACTION_MODE);

  ttkSetEnumMacro(ValidationMode, VALIDATION_MODE);
  vtkGetEnumMacro(ValidationMode, VALIDATION_MODE);

  ttkSetEnumMacro(CellMode, CELL_MODE);
  vtkGetEnumMacro(CellMode, CELL_MODE);

  vtkSetMacro(MaskOnly, bool);
  vtkGetMacro(MaskOnly, bool);

  vtkSetMacro(OutputType, int);
  vtkGetMacro(OutputType, int);

  vtkSetMacro(ExpressionString, const std::string &);
  vtkGetMacro(ExpressionString, std::string);

  vtkSetMacro(ExtractUniqueValues, bool);
  vtkGetMacro(ExtractUniqueValues, bool);

  vtkSetMacro(ArrayAttributeType, int);
  vtkGetMacro(ArrayAttributeType, int);

  vtkSetMacro(OutputArrayName, const std::string &);
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

  int AddMaskArray(vtkDataObject *output,
                   vtkDataObject *input,
                   const std::vector<double> &labels);

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
  vtkTypeMacro(ttkExtract, ttkAlgorithm);

protected:
  ttkExtract();
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
