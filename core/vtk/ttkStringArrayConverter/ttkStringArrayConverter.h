/// \class ttkStringArrayConverter
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output content of the data.csv file of the database in form of a
/// vtkTable

#pragma once

// Module include
#include <ttkStringArrayConverterModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKSTRINGARRAYCONVERTER_EXPORT ttkStringArrayConverter
  : public ttkAlgorithm {

public:
  static ttkStringArrayConverter *New();
  vtkTypeMacro(ttkStringArrayConverter, ttkAlgorithm);

protected:
  ttkStringArrayConverter();
  ~ttkStringArrayConverter() = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
