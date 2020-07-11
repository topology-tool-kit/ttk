/// \ingroup vtk
/// \class ttkOFFWriter
/// \author Julien Tierny
/// \date December 2017.<julien.tierny@lip6.fr>
/// \brief ttkOFFWriter - Object File Format Writer
///
/// Writes an .off file into VTK format.

#pragma once

#include <vtkDataSetWriter.h>

#include <Debug.h>
#include <ttkOFFWriterModule.h>

#include <fstream>

class TTKOFFWRITER_EXPORT ttkOFFWriter : public vtkDataSetWriter,
                                         protected ttk::Debug {

public:
  vtkTypeMacro(ttkOFFWriter, vtkDataSetWriter);

  static ttkOFFWriter *New();

  void PrintSelf(std::ostream &os, vtkIndent indent) override;

  // Description:
  // Specify file name of the .off file.
  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

protected:
  ttkOFFWriter();
  ~ttkOFFWriter();

  int OpenFile();
  virtual void WriteData() override;

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkOFFWriter(const ttkOFFWriter &) = delete;
  void operator=(const ttkOFFWriter &) = delete;
};
