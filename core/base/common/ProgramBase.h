/// \ingroup base
/// \class ttk::ProgramBase
/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Base editor class for standalone programs. This class parses the
/// the comamnd line, execute the TTK module and takes care of the IO.

#ifndef EDITOR_BASE_H
#define EDITOR_BASE_H

// base code includes
#include <CommandLineParser.h>

namespace ttk {

  class ProgramBase : public Debug {

  public:
    ProgramBase() {

      lastObject_ = true;
      ttk::globalDebugLevel_ = 3;

      outputPath_ = "output";
      ttkModule_ = NULL;
    }

    virtual ~ProgramBase(){};

    virtual int init(int &argc, char **argv) {

      if(!ttkModule_)
        return -1;

      std::vector<std::string> inputPaths;

      parser_.setArgument(
        "i", &inputPaths, "Input data-sets (*.vti, *vtu, *vtp)");
      parser_.setArgument(
        "o", &outputPath_, "Output file name base (no extension)", true);

      parser_.parse(argc, argv);
      debugLevel_ = ttk::globalDebugLevel_;

      threadNumber_ = ttk::globalThreadNumber_;

      int ret = 0;
      ret = load(inputPaths);

      if(ret != 0)
        return ret;

      return 0;
    }

    virtual int run() {

      if(!ttkModule_)
        return -1;

      ttkModule_->setDebugLevel(debugLevel_);
      ttkModule_->setThreadNumber(threadNumber_);

      return execute();
    }

    /// Save the output(s) of the TTK module.
    virtual int save() const = 0;

    CommandLineParser parser_;

  protected:
    std::string outputPath_;

    Debug *ttkModule_;

    /// Execute your TTK module here.
    virtual int execute() = 0;

    /// Load a sequence of input data-sets.
    virtual int load(const std::vector<std::string> &inputPaths) = 0;
  };

  template <class ttkModule>
  class Program : public ProgramBase {

  public:
    Program() {
      ttkModule_ = new ttkModule;
    }

    ~Program() {
      if(ttkModule_)
        delete ttkModule_;
    }
  };
} // namespace ttk

#endif // EDITOR_BASE_H
