/// \ingroup vtk
/// \class ttkScalarFieldNormalizer
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2016.
///
/// \brief TTK VTK-filter that normalizes an input scalar field.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output normalized scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa vtkCellDataConverter
/// \sa vtkPointDataConverter
/// \sa vtkScalarFieldSmoother
#ifndef _TTK_SCALARFIELDNORMALIZER_H
#define _TTK_SCALARFIELDNORMALIZER_H

#ifndef _MSC_VER
// ttk code includes
#include                  <Wrapper.h>

// VTK includes -- to adapt
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
#include                  <vtkSmartPointer.h>
#else
// VTK includes -- to adapt
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
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <Wrapper.h>
#endif

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkScalarFieldNormalizer
#else
class ttkScalarFieldNormalizer
#endif
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkScalarFieldNormalizer* New();
    
    vtkTypeMacro(ttkScalarFieldNormalizer, vtkDataSetAlgorithm);
    
    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
      Modified();
    }
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }   
    
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    
  protected:
    
    ttkScalarFieldNormalizer();
    
    ~ttkScalarFieldNormalizer();
    
    int normalize(vtkDataArray *input, vtkDataArray *output) const;
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
    string                ScalarField;
    vtkDataArray          *outputScalarField_;
    
    // base code features
    int doIt(vtkDataSet *input, vtkDataSet *output);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _TTK_SCALARFIELDNORMALIZER_H
