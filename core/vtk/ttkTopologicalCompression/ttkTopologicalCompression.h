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

#pragma once

// ttk code includes
#include <TopologicalCompression.h>
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkTopologicalCompressionModule.h>

class TTKTOPOLOGICALCOMPRESSION_EXPORT ttkTopologicalCompression
  : public ttkAlgorithm,
    protected ttk::TopologicalCompression {

public:
  static ttkTopologicalCompression *New();
  vtkTypeMacro(ttkTopologicalCompression, ttkAlgorithm);

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
  ttkTopologicalCompression();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string ScalarField{};
  int ScalarFieldId{0};
};
