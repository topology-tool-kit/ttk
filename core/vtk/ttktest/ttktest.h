#pragma once

// VTK Module
#include <ttktestModule.h>

// VTK Includes
#include <ttkAlgorithm.h>



// TTK Base Includes
#include <test.h>

class TTKTEST_EXPORT ttktest
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::test // and we inherit from the base class
{


public:

  static ttktest *New();
  vtkTypeMacro(ttktest, ttkAlgorithm);

protected:


  ttktest();
  ~ttktest() override = default;
 
  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation* info) override;
  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
};
