/// \ingroup vtk
/// \class ttkDistanceField
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK VTK-filter for distance field computations.
///
/// This filter takes a list of sources (a set of points with their global
/// identifiers attached to them) and produces a distance field to the closest
/// source.
///
/// \param Input0 Input geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Input1 Input sources (vtkPointSet)
/// \param Output Output distance field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "A note on two problems in connexion with graphs" \n
/// Edsger W. Dijkstra \n
/// Numerische Mathematik, 1959.
///
/// \sa ttk::DistanceField.cpp
/// \sa vtkIdentifiers
///
#pragma once

// VTK includes
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

// VTK Module
#include <ttkDistanceFieldModule.h>

// ttk code includes
#include <DistanceField.h>
#include <ttkAlgorithm.h>

enum DistanceType { Float = 0, Double = 1 };

class TTKDISTANCEFIELD_EXPORT ttkDistanceField : public ttkAlgorithm,
                                                 protected ttk::DistanceField {
public:
  static ttkDistanceField *New();

  vtkTypeMacro(ttkDistanceField, ttkAlgorithm);

  vtkSetMacro(OutputScalarFieldType, int);
  vtkGetMacro(OutputScalarFieldType, int);

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(ForceInputVertexScalarField, bool);
  vtkGetMacro(ForceInputVertexScalarField, bool);

protected:
  ttkDistanceField();
  ~ttkDistanceField() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputVertexScalarField{false};
  int OutputScalarFieldType{DistanceType::Float};
  std::string OutputScalarFieldName{"DistanceFieldValues"};
};
