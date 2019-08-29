/// \ingroup base
/// \class ttk::CommandLineParser
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2013.
///
/// \brief Basic command line parsing.

#ifndef _COMMAND_LINE_PARSER_H
#define _COMMAND_LINE_PARSER_H

#include <Debug.h>

namespace ttk {

  class CommandLineParser : public Debug {

  public:
    // 1) constructors, destructors, operators, etc.
    class CommandLineArgument : public Debug {

    public:
      CommandLineArgument() {
        boolValue_ = NULL;
        intValue_ = NULL;
        intValueList_ = NULL;
        doubleValue_ = NULL;
        doubleValueList_ = NULL;
        stringValue_ = NULL;
        stringValueList_ = NULL;
        isSet_ = false;
      };

      int print(std::stringstream &s) const {
        s << "[CommandLine]   ";
        if((isAnOption_) || (isOptional_)) {
          s << "[";
        }

        s << "-" << key_;

        if(isAnOption_)
          s << ":";

        s << " ";
        if(!isAnOption_) {
          s << "<";
          if((stringValueList_) || (intValueList_) || (doubleValueList_)) {
            s << "{";
          }
        }

        if(description_.length()) {
          s << description_;
        } else {
          s << "no description";
        }

        if((stringValue_) || (intValue_) || (doubleValue_)) {
          s << " (default: ";
          if(stringValue_)
            s << "`" << *stringValue_ << "'";
          if(intValue_)
            s << *intValue_;
          if(doubleValue_)
            s << *doubleValue_;
          s << ")";
        }

        if(isAnOption_) {
          s << " (default: " << *boolValue_ << ")";
        }

        if(!isAnOption_) {
          if((stringValueList_) || (intValueList_) || (doubleValueList_)) {
            s << "}";
          }
          s << ">";
        }

        if((isAnOption_) || (isOptional_)) {
          s << "]";
        }
        s << std::endl;

        return 0;
      };

      bool isOptional_, isAnOption_, isSet_;
      bool *boolValue_;
      int *intValue_;
      double *doubleValue_;
      std::string *stringValue_;
      std::vector<int> *intValueList_;
      std::vector<double> *doubleValueList_;
      std::vector<std::string> *stringValueList_;

      std::string key_;
      std::string description_;
    };

    CommandLineParser() {
      setArgument("d", &(ttk::globalDebugLevel_), "Global debug level", true);
      setArgument("t", &ttk::globalThreadNumber_, "Global thread number", true);
    };

    ~CommandLineParser(){};

    // 2) functions
    int parse(int argc, char **argv) {

      for(int i = 0; i < argc; i++) {

        if((std::string(argv[i]) == "-h")
           || (std::string(argv[i]) == "--help")) {
          printUsage(argv[0]);
        }

        for(int j = 0; j < (int)arguments_.size(); j++) {

          if(!arguments_[j].isAnOption_) {
            if((std::string(argv[i]) == "-" + arguments_[j].key_)
               && (i + 1 < argc)) {

              if(arguments_[j].stringValue_) {
                // let's process a std::string
                (*arguments_[j].stringValue_) = argv[i + 1];
                arguments_[j].isSet_ = true;
              } else if(arguments_[j].stringValueList_) {
                arguments_[j].stringValueList_->push_back(argv[i + 1]);
                arguments_[j].isSet_ = true;
              } else if(arguments_[j].intValue_) {
                std::stringstream s(argv[i + 1]);
                s >> *(arguments_[j].intValue_);
                arguments_[j].isSet_ = true;
              } else if(arguments_[j].intValueList_) {
                std::stringstream s(argv[i + 1]);
                arguments_[j].intValueList_->resize(
                  arguments_[j].intValueList_->size() + 1);
                s >> arguments_[j].intValueList_->back();
                arguments_[j].isSet_ = true;
              } else if(arguments_[j].doubleValue_) {
                std::stringstream s(argv[i + 1]);
                s >> *(arguments_[j].doubleValue_);
                arguments_[j].isSet_ = true;
              } else if(arguments_[j].doubleValueList_) {
                std::stringstream s(argv[i + 1]);
                arguments_[j].doubleValueList_->resize(
                  arguments_[j].doubleValueList_->size() + 1);
                s >> arguments_[j].doubleValueList_->back();
                arguments_[j].isSet_ = true;
              }
            }
          } else {
            if(std::string(argv[i]) == "-" + arguments_[j].key_) {
              *(arguments_[j].boolValue_) = !(*(arguments_[j].boolValue_));
            }
          }
        }
      }

      // check all the necessary arguments have been provided
      for(int i = 0; i < (int)arguments_.size(); i++) {
        if(!arguments_[i].isOptional_) {
          if(!arguments_[i].isSet_) {
            std::stringstream msg;
            msg << "[CommandLine] Missing mandatory argument:" << std::endl;
            arguments_[i].print(msg);
            dMsg(std::cerr, msg.str(), 1);
            printUsage(argv[0]);
          }
        }
      }

      return 0;
    };

    int printArgs(std::ostream &o = std::cout) const {

      o << "[CommandLine] Options and arguments:" << std::endl;
      for(int i = 0; i < (int)arguments_.size(); i++) {
        o << "[CommandLine]   -" << arguments_[i].key_;
        o << ": ";

        if(arguments_[i].isAnOption_) {
          if(arguments_[i].boolValue_) {
            if(*(arguments_[i].boolValue_))
              o << "true";
            else
              o << "false";
          } else {
            o << "(not set)";
          }
        } else if(arguments_[i].stringValue_) {
          if(arguments_[i].isSet_) {
            o << *(arguments_[i].stringValue_);
          } else {
            o << "(not set)";
          }
        } else if(arguments_[i].stringValueList_) {
          if(!arguments_[i].isSet_) {
            o << "(not set)";
          } else {
            for(int j = 0; j < (int)arguments_[i].stringValueList_->size();
                j++) {
              o << (*(arguments_[i].stringValueList_))[j] << " ";
            }
          }
        } else if(arguments_[i].intValue_) {
          if(!arguments_[i].isSet_) {
            o << "(not set)";
          } else {
            o << *(arguments_[i].intValue_);
          }
        } else if(arguments_[i].intValueList_) {
          if(!arguments_[i].isSet_) {
            o << "(not set)";
          } else {
            for(int j = 0; j < (int)arguments_[i].intValueList_->size(); j++) {
              o << (*(arguments_[i].intValueList_))[j] << " ";
            }
          }
        } else if(arguments_[i].doubleValue_) {
          if(!arguments_[i].isSet_) {
            o << "(not set)";
          } else {
            o << *(arguments_[i].doubleValue_);
          }
        } else if(arguments_[i].doubleValueList_) {
          if(!arguments_[i].isSet_) {
            o << "(not set)";
          } else {
            for(int j = 0; j < (int)arguments_[i].doubleValueList_->size();
                j++) {
              o << (*(arguments_[i].doubleValueList_))[j] << " ";
            }
          }
        }

        o << std::endl;
      }

      return 0;
    }

    int printUsage(const std::string &binPath) const {

      std::stringstream msg;
      msg << "[CommandLine]" << std::endl;
      msg << "[CommandLine] Usage:" << std::endl;
      msg << "[CommandLine]   " << binPath << std::endl;
      msg << "[CommandLine] Argument(s):" << std::endl;
      for(int i = 0; i < (int)arguments_.size(); i++) {
        if(!arguments_[i].isAnOption_) {
          arguments_[i].print(msg);
        }
      }
      msg << "[CommandLine] Option(s):" << std::endl;
      for(int i = 0; i < (int)arguments_.size(); i++) {
        if(arguments_[i].isAnOption_) {
          arguments_[i].print(msg);
        }
      }

      dMsg(std::cerr, msg.str(), 1);

      exit(0);
      return 0;
    };

    int setOption(const std::string &key,
                  bool *value,
                  const std::string &description = "") {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = true;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().boolValue_ = value;
      arguments_.back().isAnOption_ = true;

      return 0;
    };

    int setArgument(const std::string &key,
                    double *value,
                    const std::string &description = "",
                    const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().doubleValue_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

    int setArgument(const std::string &key,
                    std::vector<double> *value,
                    const std::string &description = "",
                    const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().doubleValueList_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

    inline int setArgument(const std::string &key,
                           int *value,
                           const std::string &description = "",
                           const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().intValue_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

    int setArgument(const std::string &key,
                    std::vector<int> *value,
                    const std::string &description = "",
                    const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().intValueList_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

    int setArgument(const std::string &key,
                    std::string *value,
                    const std::string &description = "",
                    const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().stringValue_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

    int setArgument(const std::string &key,
                    std::vector<std::string> *value,
                    const std::string &description = "",
                    const bool &optional = false) {

      if(!value)
        return -1;

      arguments_.resize(arguments_.size() + 1);
      arguments_.back().isOptional_ = optional;
      arguments_.back().key_ = key;
      arguments_.back().description_ = description;
      arguments_.back().stringValueList_ = value;
      arguments_.back().isAnOption_ = false;

      return 0;
    };

  protected:
    std::vector<CommandLineArgument> arguments_;
  };
} // namespace ttk

#endif
