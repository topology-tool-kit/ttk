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

// VTK Module
#include <ttkImportEmbeddingFromTableModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKIMPORTEMBEDDINGFROMTABLE_EXPORT ttkImportEmbeddingFromTable
  : public ttkAlgorithm {

public:
  static ttkImportEmbeddingFromTable *New();
  vtkTypeMacro(ttkImportEmbeddingFromTable, ttkAlgorithm)

    vtkSetMacro(XColumn, const std::string &);
  vtkGetMacro(XColumn, std::string);

  vtkSetMacro(YColumn, const std::string &);
  vtkGetMacro(YColumn, std::string);

  vtkSetMacro(ZColumn, const std::string &);
  vtkGetMacro(ZColumn, std::string);

  vtkSetMacro(Embedding2D, bool);
  vtkGetMacro(Embedding2D, bool);

protected:
  ttkImportEmbeddingFromTable() {
    this->setDebugMsgPrefix("ImportEmbeddingFromTable");
    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(1);
  }

  ~ttkImportEmbeddingFromTable() override{};

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  std::string XColumn;
  std::string YColumn;
  std::string ZColumn;
  bool Embedding2D;
};
