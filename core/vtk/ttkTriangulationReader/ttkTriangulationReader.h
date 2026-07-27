/// \class ttkTriangulationReader
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date April 2021
///
/// \brief TTK VTK-filter that reads a TTK Triangulation file.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output Preconditioned triangulation attached to the input dataset

#pragma once

// Module include
#include <ttkTriangulationReaderModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKTRIANGULATIONREADER_EXPORT ttkTriangulationReader
  : public ttkAlgorithm {

public:
  static ttkTriangulationReader *New();
  vtkTypeMacro(ttkTriangulationReader, ttkAlgorithm);

  vtkSetMacro(TriangulationFilePath, const std::string &);
  vtkGetMacro(TriangulationFilePath, std::string);

protected:
  ttkTriangulationReader();

  int validateFilePath();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string TriangulationFilePath{""};
};
