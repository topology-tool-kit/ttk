/// \ingroup vtk
/// \class ttkGridLayout
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.10.2018
///
/// \brief TTK VTK-filter that arranges vtkDataSets on a grid.
///
/// This filter arranges vtkDataSets which are stored in a vtkMultiBlockDataSet
/// on a grid.
///
/// \param Input vtkMultiBlockDataSet
/// \param Output vtkMultiBlockDataSet
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK Module
#include <ttkGridLayoutModule.h>

// TTK includes
#include <ttkAlgorithm.h>

class TTKGRIDLAYOUT_EXPORT ttkGridLayout : public ttkAlgorithm {

private:
  int ColAxis{1};
  int RowAxis{0};

  double ColGap{0};
  double RowGap{0};

  int NumberOfRows{0};

public:
  static ttkGridLayout *New();
  vtkTypeMacro(ttkGridLayout, ttkAlgorithm);

  vtkSetMacro(ColAxis, int);
  vtkGetMacro(ColAxis, int);

  vtkSetMacro(RowAxis, int);
  vtkGetMacro(RowAxis, int);

  vtkSetMacro(ColGap, double);
  vtkGetMacro(ColGap, double);

  vtkSetMacro(RowGap, double);
  vtkGetMacro(RowGap, double);

  vtkSetMacro(NumberOfRows, int);
  vtkGetMacro(NumberOfRows, int);

protected:
  ttkGridLayout();
  ~ttkGridLayout() override;

  int CopyObject(vtkDataObject *output, vtkDataObject *input);
  int TranslateObject(vtkDataObject *input,
                      const size_t &colAxis,
                      const size_t &rowAxis,
                      const double &dw,
                      const double &dh);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
