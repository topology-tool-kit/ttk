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

#include <ttkOFFReaderModule.h>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>

#include <string>
#include <vector>

class vtkPoints;
class vtkDoubleArray;

class TTKOFFREADER_EXPORT ttkOFFReader : public vtkUnstructuredGridAlgorithm {
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

  int countVertsData(std::string line);
  int countCellsData(std::string line);
  int processLineVert(vtkIdType curLine, std::string &line);
  int processLineCell(vtkIdType curLine, std::string &line);

private:
  ttkOFFReader(const ttkOFFReader &) = delete;
  void operator=(const ttkOFFReader &) = delete;

  char *FileName{};
  vtkIdType nbVerts_{}, nbCells_{};
  vtkIdType nbVertsData_{}, nbCellsData_{};

  vtkNew<vtkUnstructuredGrid> mesh_{};
  vtkNew<vtkPoints> points_{};
  std::vector<vtkNew<vtkDoubleArray>> vertScalars_{};
  std::vector<vtkNew<vtkDoubleArray>> cellScalars_{};
};
