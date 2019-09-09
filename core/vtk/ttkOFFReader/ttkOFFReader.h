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

#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include <vtkFiltersCoreModule.h>

#include <string>
#include <vector>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkOFFReader
#else
class ttkOFFReader
#endif
  : public vtkUnstructuredGridAlgorithm {
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

  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *) override;

  int countVertsData(std::string line);
  int countCellsData(std::string line);
  int processLineVert(vtkIdType curLine, std::string &line);
  int processLineCell(vtkIdType curLine, std::string &line);

private:
  ttkOFFReader(const ttkOFFReader &) = delete;
  void operator=(const ttkOFFReader &) = delete;

  char *FileName;
  vtkIdType nbVerts_, nbCells_;
  vtkIdType nbVertsData_, nbCellsData_;

  vtkSmartPointer<vtkUnstructuredGrid> mesh_;
  vtkSmartPointer<vtkPoints> points_;
  std::vector<vtkSmartPointer<vtkDoubleArray>> vertScalars_;
  std::vector<vtkSmartPointer<vtkDoubleArray>> cellScalars_;
};
