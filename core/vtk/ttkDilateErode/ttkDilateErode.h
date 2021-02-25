/// \ingroup vtk
/// \class ttkDilateErode
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.02.2019
///
/// \brief TTK VTK-filter that dilates or erodes a specified vertex label.
///
/// VTK wrapping code for the @DilateErode package.
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

/// \sa ttk::DilateErode

#pragma once

// VTK Module
#include <ttkDilateErodeModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <DilateErode.h>

class TTKDILATEERODE_EXPORT ttkDilateErode : public ttkAlgorithm,
                                             public ttk::DilateErode {

private:
  int Mode{0};
  std::string PivotLabel{"0"};
  int Iterations{1};

public:
  static ttkDilateErode *New();
  vtkTypeMacro(ttkDilateErode, ttkAlgorithm);

  vtkSetMacro(Mode, int);
  vtkGetMacro(Mode, int);
  vtkSetMacro(PivotLabel, std::string);
  vtkGetMacro(PivotLabel, std::string);
  vtkSetMacro(Iterations, int);
  vtkGetMacro(Iterations, int);

protected:
  ttkDilateErode();
  ~ttkDilateErode();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
