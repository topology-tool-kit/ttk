/// \ingroup vtk
/// \class ttkPointDataConverter
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK VTK-filter that converts data types for point-based scalar
/// fields (for instance, from double to float).
///
/// \param Input Input point-based scalar field (vtkDataSet)
/// \param Output Output point-based scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkCellDataConverter
///

#pragma once

// VTK Module
#include <ttkPointDataConverterModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKPOINTDATACONVERTER_EXPORT ttkPointDataConverter : public ttkAlgorithm {

public:
  static ttkPointDataConverter *New();

  vtkTypeMacro(ttkPointDataConverter, ttkAlgorithm);

  void SetOutputType(int outputType) {
    OutputType = static_cast<SupportedType>(outputType);
    Modified();
  }
  int GetOutputType() {
    return static_cast<int>(OutputType);
  }

  vtkGetMacro(UseNormalization, bool);
  vtkSetMacro(UseNormalization, bool);

protected:
  ttkPointDataConverter();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <typename InputFieldType,
            typename OutputFieldType,
            typename OutputVTKArrayType>
  int convert(vtkDataArray *inputData, vtkDataSet *output);

private:
  enum class SupportedType {
    Char = 0,
    Double,
    Float,
    Int,
    IdType,
    Short,
    UnsignedShort,
    UnsignedChar,
  };

  SupportedType OutputType{SupportedType::Char};
  bool UseNormalization{false};
};
