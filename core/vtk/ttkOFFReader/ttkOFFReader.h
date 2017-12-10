// .NAME ttkOFFReader - Object File Format Reader
// .SECTION Load a .off file into VTK format
// Note: This reader is not able to deal with comment on the file

#pragma once

#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridAlgorithm.h"

#include <string>
#include <vector>

class ttkOFFReader : public vtkUnstructuredGridAlgorithm
{
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

   int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

   int countNbFields(std::string line);
   int processLineVert(int curLine, std::string &line);
   int processLineCell(int curLine, std::string &line);

  private:
   ttkOFFReader(const ttkOFFReader &) = delete;
   void operator=(const ttkOFFReader &) = delete;

   char *FileName;
   int   nbFields_;

   vtkSmartPointer<vtkUnstructuredGrid>         mesh_;
   vtkSmartPointer<vtkPoints>                   points_;
   std::vector<vtkSmartPointer<vtkDoubleArray>> scalars_;
};
