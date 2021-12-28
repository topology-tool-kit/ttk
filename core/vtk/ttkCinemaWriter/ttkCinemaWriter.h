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
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK includes
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

// VTK Module
#include <ttkCinemaWriterModule.h>

// TTK Writer
#include <ttkTopologicalCompressionWriter.h>

class TTKCINEMAWRITER_EXPORT ttkCinemaWriter : public ttkAlgorithm {
public:
  enum class FORMAT { VTK = 0, PNG = 1, TTK = 2 };

private:
  std::string DatabasePath{""};
  int CompressionLevel{5};
  bool IterateMultiBlock{true};
  bool ForwardInput{true};
  FORMAT Format{FORMAT::VTK};

  // topological compression
  double Tolerance{1.0};
  double MaximumError{};
  double ZFPTolerance{50};
  int CompressionType{
    static_cast<int>(ttk::CompressionType::PersistenceDiagram)};
  int SQMethodPV{};
  bool ZFPOnly{false};
  bool Subdivide{false};
  bool UseTopologicalSimplification{true};

public:
  static ttkCinemaWriter *New();
  vtkTypeMacro(ttkCinemaWriter, ttkAlgorithm);

  vtkSetMacro(DatabasePath, const std::string &);
  vtkGetMacro(DatabasePath, std::string);

  ttkSetEnumMacro(Format, FORMAT);
  vtkGetEnumMacro(Format, FORMAT);

  vtkSetMacro(CompressionLevel, int);
  vtkGetMacro(CompressionLevel, int);

  vtkSetMacro(IterateMultiBlock, bool);
  vtkGetMacro(IterateMultiBlock, bool);

  vtkSetMacro(ForwardInput, bool);
  vtkGetMacro(ForwardInput, bool);

  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(MaximumError, double);
  vtkSetMacro(MaximumError, double);
  vtkGetMacro(ZFPTolerance, double);
  vtkSetMacro(ZFPTolerance, double);
  vtkGetMacro(ZFPOnly, bool);
  vtkSetMacro(ZFPOnly, bool);
  vtkGetMacro(CompressionType, int);
  vtkSetMacro(CompressionType, int);
  vtkGetMacro(Subdivide, bool);
  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);
  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkSetMacro(SQMethodPV, int);

  int DeleteDatabase();
  int GetLockFilePath(std::string &path);
  int InitializeLockFile();

protected:
  ttkCinemaWriter();
  ~ttkCinemaWriter();

  int ValidateDatabasePath();
  int ProcessDataProduct(vtkDataObject *input);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
