/// \ingroup vtk
/// \class ttkFlattenMultiBlock
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date September 2021
///
/// \brief TTK processing package for flattening the top-level
/// hierarchy of a tree vtkMultiBlockDataSet structure.

#pragma once

#include <ttkFlattenMultiBlockModule.h>

#include <ttkAlgorithm.h>

class TTKFLATTENMULTIBLOCK_EXPORT ttkFlattenMultiBlock : public ttkAlgorithm {

public:
  static ttkFlattenMultiBlock *New();

  vtkTypeMacro(ttkFlattenMultiBlock, ttkAlgorithm);

protected:
  ttkFlattenMultiBlock();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
