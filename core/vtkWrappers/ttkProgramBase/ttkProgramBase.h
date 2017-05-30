/// \ingroup vtkWrappers
/// \class ttkProgramBase
/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Base VTK editor class for standalone programs. This class parses the 
/// the comamnd line, execute the TTK module and takes care of the IO.

#ifndef _TTK_EDITOR_BASE_H
#define _TTK_EDITOR_BASE_H

// base code includes
#include                  <ProgramBase.h>

#include                  <ttkWrapper.h>

// VTK IO
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkImageData.h>
#include                  <vtkPolyData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkXMLImageDataReader.h>
#include                  <vtkXMLImageDataWriter.h>
#include                  <vtkXMLPolyDataReader.h>
#include                  <vtkXMLPolyDataWriter.h>
#include                  <vtkXMLUnstructuredGridReader.h>
#include                  <vtkXMLUnstructuredGridWriter.h>

class ttkProgramBase : public ProgramBase{

  public:
      
    ttkProgramBase(){
      
      vtkWrapper_ = NULL;
    };
      
    ~ttkProgramBase(){};
   
    /// Set the arguments of your ttk module and execute it here.
    int execute();

    vtkDataSet* getInput(const int &inputId){
      if((inputId < 0)||(inputId >= (int) inputs_.size()))
        return NULL;
      return inputs_[inputId];
    }

    int getNumberOfInputs() { return inputs_.size();};

    virtual int run(){
      
      if(!vtkWrapper_){
        return -1;
      }
        
      return execute();
    }
    
    /// Save the output(s) of the TTK module.
    virtual int save() const;
    
    virtual int setTTKmodule(vtkDataSetAlgorithm *ttkModule){
      
      vtkWrapper_ = ttkModule;
      ttkModule_ = (Debug *) ttkModule;
      
      return 0;
    }

  protected:
    
    vector<vtkDataSet *>          inputs_;
    vector<vtkSmartPointer<vtkXMLImageDataReader>>
                                  imageDataReaders_;
    vector<vtkSmartPointer<vtkXMLPolyDataReader> >
                                  polyDataReaders_;
    vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> >
                                  unstructuredGridReaders_;
    vtkDataSetAlgorithm           *vtkWrapper_;
                                  
        
    template <class vtkReaderClass>
      int load(
        const string &fileName,
        vector<vtkSmartPointer<vtkReaderClass>> &readerList);
        
    /// Load a sequence of input data-sets.
    virtual int load(const vector<string> &inputPaths);
    
    template <class vtkWriterClass>
      int save(const int &outputPortId) const;
   
};

template <class ttkModule> 
  class vtkProgram : public ttkProgramBase{
  
  public:
    
    vtkProgram(){
      ttkObject_ = vtkSmartPointer<ttkModule>::New();
      vtkWrapper_ = (vtkDataSetAlgorithm *) ttkObject_.GetPointer();
      ttkModule_ = (Debug *) ttkObject_.GetPointer();
    }
    
    virtual int run(){
      
      ttkObject_->setDebugLevel(globalDebugLevel_);
      ttkObject_->setThreadNumber(parser_.getThreadNumber());
      
      return ttkProgramBase::run();
    }
    
    vtkSmartPointer<ttkModule>    ttkObject_;
    
};

// template
#include                          <ttkProgramBase.cpp>

#endif // VTK_EDITOR_BASE_H
