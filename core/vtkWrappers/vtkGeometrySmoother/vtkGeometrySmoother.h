/// \ingroup vtkWrappers
/// \class vtkGeometrySmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter for geometry smoothing.
///
/// This filter is a dummy example for the development of TTK packages. It 
/// smooths an input mesh by average the vertex locations on the link of each 
/// vertex.
///
/// \param Input Input mesh (vtkDataSet)
/// \param Output Output mesh (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa vtkScalarFieldSmoother
/// \sa ttk::ScalarFieldSmoother

#ifndef _VTK_GEOMETRY_SMOOTHER_H
#define _VTK_GEOMETRY_SMOOTHER_H

// ttk code includes
#include                  <ScalarFieldSmoother.h>
#include                  <ttkWrapper.h>

// VTK includes
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkInformation.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkSmartPointer.h>

class VTKFILTERSCORE_EXPORT vtkGeometrySmoother 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static vtkGeometrySmoother* New();
    
    vtkTypeMacro(vtkGeometrySmoother, vtkDataSetAlgorithm);
   
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
    
    
  protected:
    
    vtkGeometrySmoother();
    
    ~vtkGeometrySmoother();
    
    TTK_SETUP();
    
    
  private:
    
    int                   NumberOfIterations;
    
    ScalarFieldSmoother   smoother_;
  
};

#endif // _VTK_GEOMETRY_SMOOTHER_H
