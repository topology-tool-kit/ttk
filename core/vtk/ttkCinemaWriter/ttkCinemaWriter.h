/// \ingroup vtk
/// \class ttkCinemaWriter
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that writes input to disk.
///
/// This filter stores the input as a VTK dataset to disk and updates the
/// data.csv file of a Cinema Spec D database.
///
/// \param Input vtkDataSet to be stored (vtkDataSet)

#pragma once

// VTK includes
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkCinemaWriterModule.h>

class TTKCINEMAWRITER_EXPORT ttkCinemaWriter : public ttkAlgorithm {

public:
  static ttkCinemaWriter *New();
  vtkTypeMacro(ttkCinemaWriter, ttkAlgorithm);

  vtkSetMacro(DatabasePath, std::string);
  vtkGetMacro(DatabasePath, std::string);

  vtkSetMacro(Mode, int);
  vtkGetMacro(Mode, int);

  vtkSetMacro(CompressionLevel, int);
  vtkGetMacro(CompressionLevel, int);

  vtkSetMacro(IterateMultiBlock, bool);
  vtkGetMacro(IterateMultiBlock, bool);

  int SetColorArray(
    int idx, int port, int connection, int fieldAssociation, const char *name) {
    this->ColorArray = {fieldAssociation, name};
    this->Modified();
    return 1;
  }

  int DeleteDatabase();

protected:
  ttkCinemaWriter();
  ~ttkCinemaWriter();

  int validateDatabasePath();
  int ProcessDataProduct(vtkDataObject *input);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);

private:
  std::string DatabasePath{""};
  int CompressionLevel{5};
  bool IterateMultiBlock{true};
  int Mode{0};
  std::pair<int, std::string> ColorArray{};
};
