/// \ingroup vtk
/// \class ttkBarycentricSubdivision
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the barycentricSubdivision processing
/// package.
///
/// VTK wrapping code for the @BarycentricSubdivision package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::BarycentricSubdivision
#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkBarycentricSubdivisionModule.h>

// TTK code includes
#include <BarycentricSubdivision.h>
#include <ttkTriangulationAlgorithm.h>

class TTKBARYCENTRICSUBDIVISION_EXPORT ttkBarycentricSubdivision
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkBarycentricSubdivision *New();
  vtkTypeMacro(ttkBarycentricSubdivision, vtkDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }
  vtkGetMacro(SubdivisionLevel, unsigned int);
  vtkSetMacro(SubdivisionLevel, unsigned int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  int FillInputPortInformation(int /*port*/, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

protected:
  ttkBarycentricSubdivision() {
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  TTK_SETUP();

  /**
   * @brief Allocate an output array of same type that input array
   */
  vtkSmartPointer<vtkDataArray>
    AllocateScalarField(vtkDataArray *const inputScalarField,
                        int ntuples) const;

  int InterpolateScalarFields(vtkUnstructuredGrid *const input,
                              vtkUnstructuredGrid *const output) const;

private:
  // number of subdivisions
  unsigned int SubdivisionLevel{1};

  // output 3D coordinates of generated points: old points first, then edge
  // middles, then triangle barycenters
  std::vector<float> points_{};
  // output triangles
  std::vector<ttk::LongSimplexId> cells_{};
  // generated point cell id
  std::vector<ttk::SimplexId> pointId_{};
  // generated points dimension: 0 vertex of parent triangulation, 1 edge
  // middle, 2 triangle barycenter
  std::vector<ttk::SimplexId> pointDim_{};

  // base worker
  ttk::BarycentricSubdivision baseWorker_{points_, cells_, pointId_, pointDim_};
};
