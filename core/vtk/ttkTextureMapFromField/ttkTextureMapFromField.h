/// \ingroup vtkWrappers
/// \class ttkTextureMapFromField
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter which generates a texture map from one or two point 
/// data scalar fields.
///
/// \param Input Input data set (vtkDataSet)
/// \param Output Output data set with texture coordinates (vtkDataSet)
///
/// This filter is useful to convert scalar fields to texture coordinates or to 
/// generate texture-based level lines out of a single scalar fields.
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa vtkProjectionFromField
///
#ifndef _TTK_TEXTURE_MAP_FROM_FIELD_H
#define _TTK_TEXTURE_MAP_FROM_FIELD_H

// ttk code includes
#include                  <Wrapper.h>

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

class VTKFILTERSCORE_EXPORT ttkTextureMapFromField 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkTextureMapFromField* New();
    
    vtkTypeMacro(ttkTextureMapFromField, vtkDataSetAlgorithm);

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
    
    vtkSetMacro(OnlyUComponent, bool);
    vtkGetMacro(OnlyUComponent, bool);

    vtkSetMacro(OnlyVComponent, bool);
    vtkGetMacro(OnlyVComponent, bool);
    
    vtkSetMacro(RepeatUTexture, bool);
    vtkGetMacro(RepeatUTexture, bool);

    vtkSetMacro(RepeatVTexture, bool);
    vtkGetMacro(RepeatVTexture, bool);
    
    vtkSetMacro(UComponent, string);
    vtkGetMacro(UComponent, string);

    vtkSetMacro(VComponent, string);
    vtkGetMacro(VComponent, string);
    
    
  protected:
    
    ttkTextureMapFromField();
    
    ~ttkTextureMapFromField();
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
    bool                  OnlyUComponent, OnlyVComponent,
                          RepeatUTexture, RepeatVTexture;
    string                UComponent, VComponent;
    vtkFloatArray         *textureCoordinates_;
    
    // base code features
    int doIt(vtkDataSet *input, vtkDataSet *output);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _TTK_TEXTURE_MAP_FROM_FIELD_H