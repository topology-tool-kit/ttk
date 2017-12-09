// .NAME ttkOFFReader - Object File Format Reader
// .SECTION Load a .off file into VTK format
// Note: This reader is not able to deal with comment on the file

#pragma once

#include "vtkUnstructuredGridAlgorithm.h"

class ttkOFFReader : public vtkUnstructuredGridAlgorithm {
public:
  vtkTypeMacro(ttkOFFReader, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  static ttkOFFReader *New();

  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  ttkOFFReader();
  ~ttkOFFReader() = default;

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

private:
  ttkOFFReader(const ttkOFFReader &) = delete;
  void operator=(const ttkOFFReader &) = delete;

  char *FileName;
};
