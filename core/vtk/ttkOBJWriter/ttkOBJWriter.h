/// \ingroup vtk
/// \class ttkOBJWriter
/// \author Julien Tierny
/// \date February 2018.<julien.tierny@lip6.fr>
/// \brief ttkOBJWriter - Object File Format Writer
///
/// Writes an .off file into VTK format.

#pragma once

#include <vtkDataSetWriter.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkOBJWriterModule.h>

#include <string>
#include <vector>

class TTKOBJWRITER_EXPORT ttkOBJWriter : public vtkDataSetWriter {

public:
  vtkTypeMacro(ttkOBJWriter, vtkDataSetWriter);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  static ttkOBJWriter *New();

  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

protected:
  ttkOBJWriter();
  ~ttkOBJWriter() override;

  int OpenFile();
  virtual void WriteData() override;

  char *Filename;
  ofstream Stream{};

private:
  ttkOBJWriter(const ttkOBJWriter &) = delete;
  void operator=(const ttkOBJWriter &) = delete;
};
