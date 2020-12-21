/// \ingroup vtk
/// \class ttkTextureMapFromField
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter which generates a texture map from one or two point
/// data scalar fields.
///
/// \param Input Input data set (vtkDataSet)
/// \param Output Output data set with texture coordinates (vtkDataSet)
///
/// This filter is useful to convert scalar fields to texture coordinates or to
/// generate texture-based level lines out of a single scalar fields.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkProjectionFromField
///

#pragma once

// VTK Module
#include <ttkTextureMapFromFieldModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKTEXTUREMAPFROMFIELD_EXPORT ttkTextureMapFromField
  : public ttkAlgorithm {

public:
  static ttkTextureMapFromField *New();

  vtkTypeMacro(ttkTextureMapFromField, ttkAlgorithm);

  vtkSetMacro(OnlyUComponent, bool);
  vtkGetMacro(OnlyUComponent, bool);

  vtkSetMacro(OnlyVComponent, bool);
  vtkGetMacro(OnlyVComponent, bool);

  vtkSetMacro(RepeatUTexture, bool);
  vtkGetMacro(RepeatUTexture, bool);

  vtkSetMacro(RepeatVTexture, bool);
  vtkGetMacro(RepeatVTexture, bool);

protected:
  ttkTextureMapFromField();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool OnlyUComponent{true}, OnlyVComponent{false}, RepeatUTexture{false},
    RepeatVTexture{false};
};
