/// \ingroup vtk
/// \class ttkCinemaImaging
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.9.2018
///
/// \brief TTK VTK-filter that generates images of a vtkDataSet.
///
/// This filter takes images of a vtkDataObject from positions specified on a
/// vtkPointSet. Each image will be a block of a vtkMultiBlockDataSet where
/// block order corresponds to point order. Each sample point can optionally
/// have vtkDoubleArrays to override the  rendering parameters, i.e, the
/// resolution, focus, clipping planes, and viewport height.
///
/// VTK wrapping code for the ttk::CinemaImaging package.
///
/// \param Input vtkDataObject that will be depicted (vtkDataObject)
/// \param Input vtkPointSet that records the camera sampling locations
/// (vtkPointSet) \param Output vtkMultiBlockDataSet that represents a list of
/// images (vtkMultiBlockDataSet)
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/geometryApproximation/">Geometry
///   Approximation example</a> \n

#pragma once

// VTK Module
#include <ttkCinemaImagingModule.h>

// TTK includes
#include <ttkAlgorithm.h>

namespace ttk {
  class CinemaImaging;
}
class vtkMultiBlockDataSet;
class vtkPointSet;
class vtkFieldData;
class vtkImageData;
class vtkPointData;
class vtkCellArray;

class TTKCINEMAIMAGING_EXPORT ttkCinemaImaging : public ttkAlgorithm {

private:
  int Backend{0};

  int Resolution[2]{256, 256};

  int ProjectionMode{0};

  bool AutoFocalPoint{true};
  bool AutoNearFar{true};
  bool AutoHeight{true};

  double FocalPoint[3]{0, 0, 0};
  double NearFar[2]{0, 1};
  double Height{1};

  double Angle{45}; // only used for perpective view

public:
  static ttkCinemaImaging *New();
  vtkTypeMacro(ttkCinemaImaging, ttkAlgorithm);

  // Backend
  vtkSetMacro(Backend, int);
  vtkGetMacro(Backend, int);

  // General Settings
  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

  // Camera
  vtkSetMacro(ProjectionMode, int);
  vtkGetMacro(ProjectionMode, int);

  vtkSetMacro(AutoFocalPoint, bool);
  vtkGetMacro(AutoFocalPoint, bool);
  vtkSetVector3Macro(FocalPoint, double);
  vtkGetVector3Macro(FocalPoint, double);

  vtkSetMacro(AutoNearFar, bool);
  vtkGetMacro(AutoNearFar, bool);
  vtkSetVector2Macro(NearFar, double);
  vtkGetVector2Macro(NearFar, double);

  vtkSetMacro(AutoHeight, bool);
  vtkGetMacro(AutoHeight, bool);
  vtkSetMacro(Height, double);
  vtkGetMacro(Height, double);

  // Perspective
  vtkSetMacro(Angle, double);
  vtkGetMacro(Angle, double);

  static vtkCellArray *GetCells(vtkPointSet *pointSet);

  static int Normalize(vtkDataArray *depthArray, const double nearFar[2]);

  static int AddFieldDataArray(vtkFieldData *fd,
                               vtkDataArray *array,
                               int tupelIdx,
                               const std::string &name = "");

  static int AddAllFieldDataArrays(vtkPointSet *inputGrid,
                                   vtkImageData *image,
                                   int tupelIdx);

  static int ComputeDirFromFocalPoint(vtkPointSet *inputGrid);

  static int EnsureGridData(vtkPointData *fd,
                            const std::string &name,
                            int nTuples,
                            const std::vector<double> &Values);

  static int MapPointAndCellData(vtkImageData *outputImage,

                                 vtkPointSet *inputObject,
                                 const ttk::CinemaImaging *renderer,
                                 const unsigned int *primitiveIdArray,
                                 const float *barycentricCoordinates,
                                 const vtkIdType *inputObjectConnectivityList);

protected:
  ttkCinemaImaging();
  ~ttkCinemaImaging();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int RequestDataSingle(vtkMultiBlockDataSet *collection,

                        vtkPointSet *object,
                        vtkPointSet *grid,
                        const std::vector<double> &defaultFocal,
                        const std::vector<double> &defaultNearFar,
                        const double defaultHeight,
                        const double defaultAngle);
};
