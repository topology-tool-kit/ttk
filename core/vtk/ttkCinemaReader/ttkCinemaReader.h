/// \class ttkCinemaReader
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

/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// Module include
#include <ttkCinemaReaderModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

class TTKCINEMAREADER_EXPORT ttkCinemaReader : public ttkAlgorithm {

public:
  static ttkCinemaReader *New();
  vtkTypeMacro(ttkCinemaReader, ttkAlgorithm);

  vtkSetMacro(DatabasePath, const std::string &);
  vtkGetMacro(DatabasePath, std::string);
  vtkSetMacro(FilePathColumnNames, const std::string &);
  vtkGetMacro(FilePathColumnNames, std::string);

protected:
  ttkCinemaReader();
  ~ttkCinemaReader() override;

  int validateDatabasePath();

  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string DatabasePath{""};
  std::string FilePathColumnNames{"FILE"};
};
