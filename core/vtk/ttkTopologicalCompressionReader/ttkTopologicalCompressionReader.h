/// \ingroup vtkWrappers
/// \class ttkTopologicalCompressionReader
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 21/04/2017
///
/// \brief VTK-filter that wraps the topologicalCompressionWriter processing
/// package.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDrivenCompression/">Persistence-Driven
///   Compression example</a> \n
///

#pragma once

// TTK
#include <TopologicalCompression.h>
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkTopologicalCompressionReaderModule.h>

class vtkImageData;

class TTKTOPOLOGICALCOMPRESSIONREADER_EXPORT ttkTopologicalCompressionReader
  : public ttkAlgorithm,
    protected ttk::TopologicalCompression {

public:
  static ttkTopologicalCompressionReader *New();

  vtkTypeMacro(ttkTopologicalCompressionReader, ttkAlgorithm);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetMacro(DataScalarType, int);
  vtkGetMacro(DataScalarType, int);

  // need this method to align with the vtkImageAlgorithm API
  vtkImageData *GetOutput();

protected:
  // Regular ImageData reader management.
  ttkTopologicalCompressionReader();
  int FillOutputPortInformation(int, vtkInformation *) override;
  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *) override;
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) override;

  // TTK management.
  void BuildMesh(vtkImageData *mesh) const;

private:
  // General properties.
  char *FileName{};

  // Data properties.
  int DataScalarType;
  std::array<int, 6> DataExtent{0, 0, 0, 0, 0, 0};
  std::array<double, 3> DataSpacing{1.0, 1.0, 1.0};
  std::array<double, 3> DataOrigin{0.0, 0.0, 0.0};
};
