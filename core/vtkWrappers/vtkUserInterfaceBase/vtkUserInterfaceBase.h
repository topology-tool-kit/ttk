/// \ingroup vtkWrappers
/// \class vtkUserInterfaceBase
/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date March 2017.
///
/// \brief Interactions and rendering.

#ifndef _VTK_USERINTERFACE_BASE_H
#define _VTK_USERINTERFACE_BASE_H

#include                  <vtkProgramBase.h>

// VTK includes
#include                  <vtkDataSetSurfaceFilter.h>
#include                  <vtkInteractorStyleSwitch.h>
#include                  <vtkInteractorStyleTrackballCamera.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPNGReader.h>
#include                  <vtkPolyDataMapper.h>
#include                  <vtkProperty.h>
#include                  <vtkRenderer.h>
#include                  <vtkRenderWindow.h>
#include                  <vtkRenderWindowInteractor.h>
#include                  <vtkTextureMapFromField.h>
#include                  <vtkVRMLExporter.h>

// forward declaration
class vtkUserInterfaceBase;

// Custom interactors
class ttkCustomInteractor : public vtkInteractorStyleTrackballCamera{

  public:
    
    static ttkCustomInteractor* New();
    
    vtkTypeMacro(ttkCustomInteractor, vtkInteractorStyleTrackballCamera);
    
    virtual void OnKeyPress();
    
    inline int setUserInterface(vtkUserInterfaceBase *userInterface){
      
      userInterface_ = userInterface;
      return 0;
    };
    
  protected:
    vtkUserInterfaceBase  *userInterface_;  
};

class ttkKeyHandler : public Debug{
  
  public:
    
    virtual int OnKeyPress(
      vtkRenderWindowInteractor *interactor, string &key) = 0;
};

class vtkUserInterfaceBase : public vtkProgramBase{
  
  public:
    
    vtkUserInterfaceBase();
    
    ~vtkUserInterfaceBase();
   
    int exportScene(const string &fileName = "output.wrl") const;
    
    ttkKeyHandler* getKeyHandler() { return keyHandler_;};

    int hideOutputs(const vector<int> &outputList){

      for(int i = 0; i < (int) outputList.size(); i++){
        bool isAlreadyIn = false;
        for(int j = 0; j < (int) hiddenOutputs_.size(); j++){
          if(outputList[i] == hiddenOutputs_[j]){
            isAlreadyIn = true;
            break;
          }
        }
        if(!isAlreadyIn)
          hiddenOutputs_.push_back(outputList[i]);
      }

      return 0;
    }

    int init(int &argc, char **argv);
  
    int refresh();

    int run();
    
    int setKeyHandler(ttkKeyHandler *handler){
      keyHandler_ = handler;
      return 0;
    }

    int switchOutput(const int &outputId);

    int switchTransparency();

  protected:
    
    bool                          hasTexture_, isUp_, repeat_, transparency_;
    vector<bool>                  visibleOutputs_;
    vector<int>                   hiddenOutputs_;
    ttkKeyHandler                 *keyHandler_;
    vector<vtkPolyData *>         surfaces_;
    vtkSmartPointer<ttkCustomInteractor>
                                  customInteractor_;
    vector<vtkSmartPointer<vtkActor> >
                                  mainActors_;
    vector<vtkSmartPointer<vtkDataSetSurfaceFilter> >
                                  boundaryFilters_;
    vtkSmartPointer<vtkPNGReader> pngReader_;
    vector<vtkSmartPointer<vtkPolyDataMapper> >
                                  boundaryMappers_;
    vtkSmartPointer<vtkRenderer>  renderer_;
    vtkSmartPointer<vtkRenderWindow>
                                  renderWindow_;
    vtkSmartPointer<vtkRenderWindowInteractor>
                                  interactor_;
    vtkSmartPointer<vtkTexture>   texture_;
    vector<vtkSmartPointer<vtkTextureMapFromField> >
                                  textureMapFromFields_;
    
    int updateScalarFieldTexture();
};

template <class ttkModule>
  class vtkUserInterface : public vtkUserInterfaceBase{
  
  public:
    
    vtkUserInterface(){
      ttkObject_ = vtkSmartPointer<ttkModule>::New();
      vtkWrapper_ = (vtkDataSetAlgorithm *) ttkObject_.GetPointer();
      ttkModule_ = (Debug *) ttkObject_.GetPointer();
    };
    
    virtual int run(){
      
      ttkObject_->setDebugLevel(globalDebugLevel_);
      ttkObject_->setThreadNumber(parser_.getThreadNumber());
      
      return vtkUserInterfaceBase::run();
    }
    
    vtkSmartPointer<ttkModule>    ttkObject_;
};


#include <vtkUserInterfaceBase.cpp>

#endif //_VTK_USERINTERFACE_BASE_H
