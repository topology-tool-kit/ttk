/// \ingroup vtkWrappers
/// \class ttkScalarFieldSmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter for scalar field smoothing.
///
/// This class is a dummy example for the development of TTK filters. It 
/// smooths an input scalar field by averaging the scalar values on the link
/// of each vertex.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
/// 
/// \sa vtkGeometrySmoother
/// \sa ttk::ScalarFieldSmoother
#ifndef _TTK_SCALAR_FIELD_SMOOTHER_H
#define _TTK_SCALAR_FIELD_SMOOTHER_H

// ttk code includes
#include                  <ScalarFieldSmoother.h>
#include                  <ttkWrapper.h>

// VTK includes
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkUnsignedShortArray.h>

class VTKFILTERSCORE_EXPORT ttkScalarFieldSmoother 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkScalarFieldSmoother* New();
    
    vtkTypeMacro(ttkScalarFieldSmoother, vtkDataSetAlgorithm);
   
    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }   
    
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
    vtkSetMacro(NumberOfIterations, int);
    vtkGetMacro(NumberOfIterations, int);
    
    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    vtkSetMacro(ScalarFieldIdentifier, int);
    vtkGetMacro(ScalarFieldIdentifier, int);

    
  protected:
    
    ttkScalarFieldSmoother();
    
    ~ttkScalarFieldSmoother();
 
    TTK_SETUP();
    
    
  private:
    
    int                   NumberOfIterations;
    int                   ScalarFieldIdentifier;
    string                ScalarField;
    vtkDataArray          *outputScalarField_;
    ScalarFieldSmoother   smoother_;
    
};

#endif // _TTK_SCALAR_FIELD_SMOOTHER_H