/// \ingroup vtk
/// \class ttkOBJWriter
/// \author Julien Tierny
/// \date February 2018.<julien.tierny@lip6.fr>
/// \brief ttkOBJWriter - Object File Format Writer
///
/// Writes an .off file into VTK format.

#pragma once

#include <vtkDataSetWriter.h>

// VTK Module
#include <Debug.h>
#include <ttkOBJWriterModule.h>

#include <fstream>

class TTKOBJWRITER_EXPORT ttkOBJWriter : public vtkDataSetWriter,
                                         protected ttk::Debug {

public:
  vtkTypeMacro(ttkOBJWriter, vtkDataSetWriter);

  static ttkOBJWriter *New();

  void PrintSelf(std::ostream &os, vtkIndent indent) override;

  // Description:
  // Specify file name of the .obj file.
  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

protected:
  ttkOBJWriter();
  ~ttkOBJWriter() override;

  int OpenFile();
  virtual void WriteData() override;

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkOBJWriter(const ttkOBJWriter &) = delete;
  void operator=(const ttkOBJWriter &) = delete;
};
