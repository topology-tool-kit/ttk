/// \ingroup vtk
/// \class ttkTriangulationWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkTriangulationWriter - Dipha Image Data Format Writer
///
/// Writes a Dipha Cubical Complex file from a VTK Image Data or a
/// Dipha Explicit Complex from a VTK Unstructured Grid dataset

#pragma once

#include <ttkAlgorithm.h>
#include <ttkTriangulationWriterModule.h>

#include <fstream>

class TTKTRIANGULATIONWRITER_EXPORT ttkTriangulationWriter
  : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkTriangulationWriter, ttkAlgorithm);

  static ttkTriangulationWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkTriangulationWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeImageData(vtkDataObject *);
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkTriangulationWriter(const ttkTriangulationWriter &) = delete;
  void operator=(const ttkTriangulationWriter &) = delete;
};
