/// \ingroup vtk
/// \class ttkImportEmbeddingFromTable
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2018
///
/// \brief TTK VTK-filter that embeds a vtkPointSet with the data of a vtkTable
///
/// \param Input Input point set (vtkPointSet)
/// \param Input Input table (vtkTable)
/// \param Output Output point set (vtkPointSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

#include <ttkImportEmbeddingFromTableModule.h>
#include <ttkAlgorithm.h>

class TTKIMPORTEMBEDDINGFROMTABLE_EXPORT ttkImportEmbeddingFromTable
  : public ttkAlgorithm {

private:
  bool Embedding2D;

public:
  vtkSetMacro(Embedding2D, bool);
  vtkGetMacro(Embedding2D, bool);

  static ttkImportEmbeddingFromTable *New();
  vtkTypeMacro(ttkImportEmbeddingFromTable, ttkAlgorithm);

protected:
  ttkImportEmbeddingFromTable();
  ~ttkImportEmbeddingFromTable();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
