/// \ingroup vtk
/// \class ttkTriangulationWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date April 2021
/// \brief ttkTriangulationWriter - Explicit Triangulation Writer
///
/// Writes the internal state of an Explicit Triangulation to disk to
/// skip preconditioning when loaded with \sa ttkTriangulationReader

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

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkTriangulationWriter(const ttkTriangulationWriter &) = delete;
  void operator=(const ttkTriangulationWriter &) = delete;
};
