/// \ingroup baseCode
/// \class ttk::ProgramBase
/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Base editor class for standalone programs. This class parses the 
/// the comamnd line, execute the TTK module and takes care of the IO.

#ifndef EDITOR_BASE_H
#define EDITOR_BASE_H

// base code includes
#include                  <CommandLineParser.h>

namespace ttk{

  class ProgramBase : public Debug{

    public:
        
      ProgramBase(){
        
        lastObject_ = true;
        globalDebugLevel_ = debugLevel_;
        
        outputPath_ = "output";
        ttkModule_ = NULL;
      }
        
      ~ProgramBase(){};
      
      virtual int init(int &argc, char **argv){
        
        if(!ttkModule_)
          return -1;
        
        vector<string> inputPaths;

        parser_.setArgument("i", &inputPaths,
          "Input data-sets (*.vti, *vtu, *vtp)");
        parser_.setArgument("o", &outputPath_,
          "Output file name base (no extension)", true);
        
        parser_.parse(argc, argv);
        debugLevel_ = globalDebugLevel_;

        threadNumber_ = parser_.getThreadNumber();
        
        int ret = 0;
        ret = load(inputPaths);
        
        if(ret != 0)
          return ret;
        
        return 0;
      }
     
      virtual int run(){
        
        if(!ttkModule_)
          return -1;
        
        ttkModule_->setDebugLevel(globalDebugLevel_);
        ttkModule_->setThreadNumber(threadNumber_);
          
        return execute();
      }
     
      /// Save the output(s) of the TTK module.
      virtual int save() const = 0;
      
      virtual int setTTKmodule(Debug *ttkModule){
        ttkModule_ = ttkModule;
        return -1;
      }
     
      CommandLineParser             parser_;
      
      
    protected:
      
      string                        outputPath_;
      
      Debug                         *ttkModule_;
     
      /// Execute your TTK module here.
      virtual int execute() = 0;
      
      /// Load a sequence of input data-sets.
      virtual int load(const vector<string> &inputPaths) = 0;
      
  };
  
  template <class ttkModule> 
    class Program : public ProgramBase{
    
    public:
      
      Program(){
        ttkModule_ = new ttkModule;
      }
      
      ~Program(){
        if(ttkModule_)
          delete ttkModule_;
      }
  };
}

#endif // EDITOR_BASE_H
