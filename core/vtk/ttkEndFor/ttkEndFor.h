/// \ingroup vtk
/// \class ttkEndFor
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that requests data as long as it is available.
///
/// This filter requests more data as long as the maximum number of elements is
/// not reached. This filter works in conjunction with the ttkForEachRow filter.
///
/// \param Input vtkDataObject that will be passed through after all iterations.
/// \param Output vtkDataObject Shallow copy of the input
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK Module
#include <ttkEndForModule.h>

// TTK includes
#include <ttkAlgorithm.h>

class TTKENDFOR_EXPORT ttkEndFor : public ttkAlgorithm {

private:
  int LastIterationIdx{-1};

public:
  static ttkEndFor *New();
  vtkTypeMacro(ttkEndFor, ttkAlgorithm);

protected:
  ttkEndFor();
  ~ttkEndFor();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};