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
#ifndef _TTK_DISTANCEFIELD_H
#define _TTK_DISTANCEFIELD_H

// VTK includes
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <DistanceField.h>
#include <ttkWrapper.h>

enum DistanceType { Float = 0, Double };

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDistanceField
#else
class ttkDistanceField
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkDistanceField *New();

  vtkTypeMacro(ttkDistanceField, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(OutputScalarFieldType, int);
  vtkGetMacro(OutputScalarFieldType, int);

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(ForceInputVertexScalarField, int);
  vtkGetMacro(ForceInputVertexScalarField, int);

  vtkSetMacro(InputVertexScalarFieldName, std::string);
  vtkGetMacro(InputVertexScalarFieldName, std::string);

  int getTriangulation(vtkDataSet *input);
  int getIdentifiers(vtkDataSet *input);

protected:
  ttkDistanceField();
  ~ttkDistanceField();

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  std::string ScalarField;
  int OutputScalarFieldType;
  std::string OutputScalarFieldName;
  bool ForceInputVertexScalarField;
  std::string InputVertexScalarFieldName;

  ttk::DistanceField distanceField_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *identifiers_;
};

#endif // _TTK_DISTANCEFIELD_H
