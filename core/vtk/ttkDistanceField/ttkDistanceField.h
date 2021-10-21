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

// VTK Module
#include <ttkDistanceFieldModule.h>

// ttk code includes
#include <DistanceField.h>
#include <ttkAlgorithm.h>

#include <string>

class TTKDISTANCEFIELD_EXPORT ttkDistanceField : public ttkAlgorithm,
                                                 protected ttk::DistanceField {
public:
  static ttkDistanceField *New();

  vtkTypeMacro(ttkDistanceField, ttkAlgorithm);

  vtkSetMacro(OutputScalarFieldType, int);
  vtkGetMacro(OutputScalarFieldType, int);

  vtkSetMacro(OutputScalarFieldName, const std::string &);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(ForceInputVertexScalarField, bool);
  vtkGetMacro(ForceInputVertexScalarField, bool);

protected:
  ttkDistanceField();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  enum class DistanceType { Float = 0, Double = 1 };

  bool ForceInputVertexScalarField{false};
  int OutputScalarFieldType{static_cast<int>(DistanceType::Float)};
  std::string OutputScalarFieldName{"DistanceFieldValues"};
};
