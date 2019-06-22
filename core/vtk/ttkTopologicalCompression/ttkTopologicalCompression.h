/// \ingroup vtkWrappers
/// \class vtkTopologicalCompression
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 20/04/2017
///
/// \brief TTK VTK-filter that wraps the topologicalCompression processing
/// package.
///
/// VTK wrapping code for the @TopologicalCompression package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::TopologicalCompression
#ifndef _VTK_TOPOLOGICALCOMPRESSION_H
#define _VTK_TOPOLOGICALCOMPRESSION_H

// ttk code includes
#include <TopologicalCompression.h>
#include <ttkWrapper.h>

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
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLImageDataWriter.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTopologicalCompression
#else
class ttkTopologicalCompression
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkTopologicalCompression *New();

  vtkTypeMacro(ttkTopologicalCompression, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  // Set/Get macros (arguments)
  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(MaximumError, double);
  vtkGetMacro(MaximumError, double);

  vtkSetMacro(CompressionType, int);
  vtkGetMacro(CompressionType, int);

  vtkSetMacro(SQMethod, std::string);
  vtkGetMacro(SQMethod, std::string);

  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(Subdivide, bool);

  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);

  inline void SetSQMethodPV(int c) {
    switch(c) {
      case 1:
        SetSQMethod("r");
        break;
      case 2:
        SetSQMethod("d");
        break;
      case 0:
      default:
        SetSQMethod("");
        break;
    }
  }

protected:
  ttkTopologicalCompression() {
    CompressionType = 0;
    Tolerance = 10;
    Subdivide = false;
    MaximumError = 10;
    UseTopologicalSimplification = true;
    ScalarFieldId = 0;
    outputScalarField_ = nullptr;
    outputOffsetField_ = nullptr;
    UseAllCores = true;
  }

  ~ttkTopologicalCompression(){};

  TTK_SETUP();

private:
  double Tolerance;
  double MaximumError;
  int CompressionType;
  std::string SQMethod;
  bool Subdivide;
  bool UseTopologicalSimplification;
  std::string ScalarField;
  int ScalarFieldId;

  vtkSmartPointer<vtkDataArray> outputScalarField_;
  vtkSmartPointer<vtkIntArray> outputOffsetField_;
  ttkTriangulation triangulation_;
  ttk::Triangulation *internalTriangulation_;
  ttk::TopologicalCompression topologicalCompression_;
};

#endif // _VTK_TOPOLOGICALCOMPRESSION_H
