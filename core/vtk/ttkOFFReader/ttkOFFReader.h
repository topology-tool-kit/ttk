/// \ingroup vtk
/// \class ttkOFFReader
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date December 2017.
/// \brief ttkOFFReader - Object File Format Reader
///
/// Load an .off file into VTK format
///
/// Note: This reader is not able to deal with comment on the file

#pragma once

#include <Debug.h>
#include <ttkOFFReaderModule.h>

#include <vtkUnstructuredGridAlgorithm.h>

class TTKOFFREADER_EXPORT ttkOFFReader : public vtkUnstructuredGridAlgorithm,
                                         protected ttk::Debug {
public:
  vtkTypeMacro(ttkOFFReader, vtkUnstructuredGridAlgorithm);

  static ttkOFFReader *New();

  void PrintSelf(std::ostream &os, vtkIndent indent) override;

  // Description:
  // Specify file name of the .off file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  ttkOFFReader();
  ~ttkOFFReader() override = default;

  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *) override;

private:
  ttkOFFReader(const ttkOFFReader &) = delete;
  void operator=(const ttkOFFReader &) = delete;

  char *FileName{};
};
