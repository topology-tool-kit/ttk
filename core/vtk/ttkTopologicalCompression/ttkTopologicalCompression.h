/// \ingroup vtkWrappers
/// \class vtkTopologicalCompression
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 20/04/2017
///
/// \brief TTK VTK-filter that wraps the topologicalCompression processing
/// package.
///
/// VTK wrapping code for the ttk::TopologicalCompression package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::TopologicalCompression
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDrivenCompression/">Persistence-Driven
///   Compression example</a> \n
///

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

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(MaximumError, double);
  vtkGetMacro(MaximumError, double);

  vtkSetMacro(CompressionType, int);
  vtkGetMacro(CompressionType, int);

  vtkSetMacro(SQMethod, const std::string &);
  vtkGetMacro(SQMethod, std::string);

  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(Subdivide, bool);

  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);

  inline void SetSQMethodPV(int c) {
    if(c == 1) {
      SetSQMethod("r");
    } else if(c == 2) {
      SetSQMethod("d");
    } else {
      SetSQMethod("");
    }
  }

protected:
  ttkTopologicalCompression();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
