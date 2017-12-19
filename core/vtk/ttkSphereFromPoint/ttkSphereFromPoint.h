/// \ingroup vtk
/// \class ttkSphereFromPoint
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter that produces sphere-only glyphs.
///
/// \param Input Input point cloud (vtkDataSet)
/// \param Output Output spheres (vtkPolyData)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
#ifndef _TTK_SPHERE_FROM_POINT_H
#define _TTK_SPHERE_FROM_POINT_H

#ifndef _MSC_VER
// ttk code includes
#include                  <Wrapper.h>

// VTK includes 
#include                  <vtkAppendPolyData.h>  
#include                  <vtkCharArray.h>
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
#include                  <vtkSphereSource.h>
#include                  <vtkType.h>
#else
// VTK includes 
#include                  <vtkAppendPolyData.h>  
#include                  <vtkCharArray.h>
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
#include                  <vtkSphereSource.h>
#include                  <vtkType.h>

// ttk code includes
#include                  <Wrapper.h>
#endif

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkSphereFromPoint
#else
class ttkSphereFromPoint
#endif 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkSphereFromPoint* New();
   
    // macros
    vtkTypeMacro(ttkSphereFromPoint, vtkDataSetAlgorithm);

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
    
    vtkSetMacro(EndPhi, int);
    
    vtkSetMacro(EndTheta, int);
    
    vtkSetMacro(PhiResolution, int);
    
    vtkSetMacro(Radius, double);
    
    vtkSetMacro(StartPhi, int);
    
    vtkSetMacro(StartTheta, int);
    
    vtkSetMacro(ThetaResolution, int);
    
    /// Over-ride the input data type to vtkDataSet.
    int FillOutputPortInformation(int port,
      vtkInformation *info){
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      return 1;
    }
    
  protected:
    
    ttkSphereFromPoint();
    
    ~ttkSphereFromPoint();
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
    int                   ThetaResolution, StartTheta, EndTheta, 
                          PhiResolution, StartPhi, EndPhi;
    double                Radius;
   
    vtkAppendPolyData     *masterAppender_;
    vector<vtkAppendPolyData *>     
                          appenderList_;
    vector<vtkSphereSource *>
                          sphereList_;
    vector<vector<vtkDataArray *> >
                          dataArrayList_;
    
    // base code features
    int doIt(vtkDataSet *input, vtkPolyData *output);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _TTK_SPHERE_FROM_POINT_H
