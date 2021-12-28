/// \ingroup vtk
/// \class ttkMorphologicalOperators
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.02.2019
///
/// \brief TTK VTK-filter that dilates or erodes a specified vertex label.
///
/// VTK wrapping code for the ttk::MorphologicalOperators package.
///
/// This filter either a) dilates a specified label by assigning the label of a
/// corresponding vertex to all its neighbors, or b) erodes a specified label by
/// assigning to a corresponding vertex the largest label among its neighbors.
///
/// The input data array that will be dilated or eroded needs to be specified
/// via the standard VTK call SetInputArrayToProcess() with the following
/// parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the VTK array)

/// \sa ttk::MorphologicalOperators

#pragma once

// VTK Module
#include <ttkMorphologicalOperatorsModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <MorphologicalOperators.h>

class TTKMORPHOLOGICALOPERATORS_EXPORT ttkMorphologicalOperators
  : public ttkAlgorithm,
    protected ttk::MorphologicalOperators {

private:
  int Mode{0};
  std::string PivotLabel{"0"};
  int Iterations{1};
  bool Grayscale{false};

public:
  static ttkMorphologicalOperators *New();
  vtkTypeMacro(ttkMorphologicalOperators, ttkAlgorithm);

  vtkSetMacro(Mode, int);
  vtkGetMacro(Mode, int);
  vtkSetMacro(PivotLabel, const std::string &);
  vtkGetMacro(PivotLabel, std::string);
  vtkSetMacro(Iterations, int);
  vtkGetMacro(Iterations, int);
  vtkSetMacro(Grayscale, bool);
  vtkGetMacro(Grayscale, bool);

protected:
  ttkMorphologicalOperators();
  ~ttkMorphologicalOperators();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
