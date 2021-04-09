/// \class ttkTriangulationReader
/// \ingroup vtk
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output content of the data.csv file of the database in form of a
/// vtkTable

#pragma once

// Module include
#include <ttkTriangulationReaderModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

class TTKTRIANGULATIONREADER_EXPORT ttkTriangulationReader
  : public ttkAlgorithm {

public:
  static ttkTriangulationReader *New();
  vtkTypeMacro(ttkTriangulationReader, ttkAlgorithm);

  vtkSetMacro(DatabasePath, std::string);
  vtkGetMacro(DatabasePath, std::string);
  vtkSetMacro(FilePathColumnNames, std::string);
  vtkGetMacro(FilePathColumnNames, std::string);

protected:
  ttkTriangulationReader();
  ~ttkTriangulationReader() override;

  int validateDatabasePath();

  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string DatabasePath{""};
  std::string FilePathColumnNames{"FILE"};
};
