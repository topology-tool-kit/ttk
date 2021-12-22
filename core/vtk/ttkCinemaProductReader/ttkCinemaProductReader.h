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
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK Module
#include <ttkCinemaProductReaderModule.h>

// VTK includes
#include <ttkAlgorithm.h>

#include <ttkTopologicalCompressionReader.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTIFFReader.h>
#include <vtkXMLGenericDataObjectReader.h>

class TTKCINEMAPRODUCTREADER_EXPORT ttkCinemaProductReader
  : public ttkAlgorithm {

public:
  static ttkCinemaProductReader *New();
  vtkTypeMacro(ttkCinemaProductReader, ttkAlgorithm);

  vtkSetMacro(FilepathColumnName, const std::string &);
  vtkGetMacro(FilepathColumnName, std::string);
  vtkSetMacro(AddFieldDataRecursively, bool);
  vtkGetMacro(AddFieldDataRecursively, bool);

protected:
  ttkCinemaProductReader();
  ~ttkCinemaProductReader();

  vtkSmartPointer<vtkDataObject> readFileLocal(const std::string &pathToFile);
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
  vtkNew<ttkTopologicalCompressionReader> topologicalCompressionReader{};

  // TIFF READER
  vtkNew<vtkTIFFReader> tiffReader{};

  // LOCAL-LEGACY && REMOTE-LEGACY
  vtkNew<vtkGenericDataObjectReader> genericDataObjectReader{};

  // LOCAL-XML
  vtkNew<vtkXMLGenericDataObjectReader> xmlGenericDataObjectReader{};
};
