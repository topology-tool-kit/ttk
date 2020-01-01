/// \ingroup vtk
/// \class ttkCinemaProductReader
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads the data products that are referenced in a
/// vtkTable.
///
/// This filter reads the products that are referenced in a vtkTable. The
/// results are stored in a vtkMultiBlockDataSet where each block corresponds to
/// a row of the table with consistent ordering.
///
/// \param Input vtkTable that contains data product references (vtkTable)
/// \param Output vtkMultiBlockDataSet where each block is a referenced product
/// of an input table row (vtkMultiBlockDataSet)

#pragma once

// VTK Module
#include <ttkCinemaProductReaderModule.h>

// VTK includes
#include <ttkAlgorithm.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkSmartPointer.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <ttkTopologicalCompressionReader.h>

class TTKCINEMAPRODUCTREADER_EXPORT ttkCinemaProductReader
  : public ttkAlgorithm {

public:
  static ttkCinemaProductReader *New();
  vtkTypeMacro(ttkCinemaProductReader, ttkAlgorithm);

  vtkSetMacro(FilepathColumnName, std::string);
  vtkGetMacro(FilepathColumnName, std::string);
  vtkSetMacro(AddFieldDataRecursively, bool);
  vtkGetMacro(AddFieldDataRecursively, bool);

protected:
  ttkCinemaProductReader();
  ~ttkCinemaProductReader();

  vtkSmartPointer<vtkDataObject> readFileLocal(std::string pathToFile);
  int addFieldDataRecursively(vtkDataObject *object, vtkFieldData *fd);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string FilepathColumnName{"FILE"};
  bool AddFieldDataRecursively{true};

  // TTK READER
  vtkSmartPointer<ttkTopologicalCompressionReader> ttkTopologicalCompressionReader_
    = vtkSmartPointer<ttkTopologicalCompressionReader>::New();

  // LOCAL-LEGACY && REMOTE-LEGACY
  vtkSmartPointer<vtkGenericDataObjectReader> genericDataObjectReader
    = vtkSmartPointer<vtkGenericDataObjectReader>::New();

  // LOCAL-XML
  vtkSmartPointer<vtkXMLGenericDataObjectReader> xmlGenericDataObjectReader
    = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
};
