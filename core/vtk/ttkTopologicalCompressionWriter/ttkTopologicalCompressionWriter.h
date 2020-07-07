/// \ingroup vtkWrappers
/// \class ttkTopologicalCompressionWriter
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 21/04/2017
///
/// \brief VTK-filter that wraps the topologicalCompressionWriter processing
/// package.

#pragma once

// TTK
#include <TopologicalCompression.h>
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkTopologicalCompressionWriterModule.h>

class vtkImageData;

class TTKTOPOLOGICALCOMPRESSIONWRITER_EXPORT ttkTopologicalCompressionWriter
  : public ttkAlgorithm,
    protected ttk::TopologicalCompression {

public:
  static ttkTopologicalCompressionWriter *New();

  vtkTypeMacro(ttkTopologicalCompressionWriter, ttkAlgorithm);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);

  vtkGetMacro(MaximumError, double);
  vtkSetMacro(MaximumError, double);

  vtkGetMacro(ZFPBitBudget, double);
  vtkSetMacro(ZFPBitBudget, double);

  vtkGetMacro(ZFPOnly, bool);
  vtkSetMacro(ZFPOnly, bool);

  vtkGetMacro(CompressionType, int);
  vtkSetMacro(CompressionType, int);

  vtkGetMacro(NbSegments, int);
  vtkSetMacro(NbSegments, int);

  vtkGetMacro(NbVertices, int);
  vtkSetMacro(NbVertices, int);

  vtkGetMacro(SQMethod, std::string);
  vtkSetMacro(SQMethod, std::string);

  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(Subdivide, bool);

  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);

  inline void SetSQMethodPV(int c) {
    if(c == 1) {
      SetSQMethod("r");
    } else if(c == 2) {
      SetSQMethod("d");
    } else if(c == 0) {
      SetSQMethod("");
    }
  }

  // expose Write method for vtkWriter API compatibility (duck-typing)
  int Write();

protected:
  // Regular writer management.
  ttkTopologicalCompressionWriter();
  virtual int FillInputPortInformation(int port, vtkInformation *info) override;

  // need those two methods to align with the vtkWriter API
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

  template <typename triangulationType>
  void PerformCompression(vtkDataArray *inputScalarField,
                          vtkDataArray *outputScalarField,
                          const triangulationType &triangulation);

private:
  // Writer parameters.
  char *FileName{};

  // Whatever.
  ttkTopologicalCompressionWriter(const ttkTopologicalCompressionWriter &);
  void operator=(const ttkTopologicalCompressionWriter &);

  // give ttkCinemaWriter access to ttkTopologicalCompressionWriter
  // protected member functions
  friend class ttkCinemaWriter;
};
