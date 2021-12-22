/// \ingroup vtk
/// \class ttkForEach
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.01.2019
///
/// \brief TTK VTK-filter that works in conjunction with the ttkEndFor filter to
/// iterate over blocks, rows, array values, and arrays.
///
/// This filter works in conjunction with the ttkEndFor filter to iterate over
/// blocks, rows, array values, and arrays.
///
/// \param Input vktObject either a vtkMultiBlockDataSet, vtkTable, or
/// vtkDataSet \param Output vktObject one element of the input
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK Module
#include <ttkForEachModule.h>

// TTK includes
#include <ttkExtract.h>

class TTKFOREACH_EXPORT ttkForEach : public ttkExtract {

private:
  vtkDataObject *LastInput{nullptr};
  int IterationIdx{0};
  int IterationNumber{0};

public:
  static ttkForEach *New();
  vtkTypeMacro(ttkForEach, ttkExtract);

  vtkGetMacro(IterationIdx, int);
  vtkSetMacro(IterationIdx, int);
  vtkGetMacro(IterationNumber, int);
  vtkSetMacro(IterationNumber, int);

protected:
  ttkForEach();
  ~ttkForEach();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};