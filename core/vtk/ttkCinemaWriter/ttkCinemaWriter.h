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
#include <vtkNew.h>

// VTK Module
#include <ttkCinemaWriterModule.h>

// TTK Writer
#include <ttkTopologicalCompressionWriter.h>

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

  vtkSetMacro(ForwardInput, bool);
  vtkGetMacro(ForwardInput, bool);

  int DeleteDatabase();

  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(MaximumError, double);
  vtkSetMacro(MaximumError, double);
  vtkGetMacro(ZFPBitBudget, double);
  vtkSetMacro(ZFPBitBudget, double);
  vtkGetMacro(ZFPOnly, bool);
  vtkSetMacro(ZFPOnly, bool);
  vtkGetMacro(CompressionType, int);
  vtkSetMacro(CompressionType, int);
  vtkGetMacro(Subdivide, bool);
  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);
  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkSetMacro(SQMethodPV, int);

protected:
  ttkCinemaWriter();
  ~ttkCinemaWriter();

  int validateDatabasePath();
  int ProcessDataProduct(vtkDataObject *input);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string DatabasePath{""};
  int CompressionLevel{5};
  bool IterateMultiBlock{true};
  bool ForwardInput{true};
  int Mode{0};

  // topological compression
  double Tolerance{1.0};
  double MaximumError{};
  double ZFPBitBudget{0};
  int CompressionType{
    static_cast<int>(ttk::CompressionType::PersistenceDiagram)};
  int SQMethodPV{};
  bool ZFPOnly{false};
  bool Subdivide{false};
  bool UseTopologicalSimplification{true};
};
