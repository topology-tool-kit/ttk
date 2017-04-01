/// \ingroup baseCode
/// \class ttk::CommandLineParser
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2013.
/// 
/// \brief Basic command line parsing.

#ifndef                 _COMMAND_LINE_PARSER_H
#define                 _COMMAND_LINE_PARSER_H

#include                <Debug.h>

namespace ttk{

  class CommandLineParser : public Debug{
    
    public:
      
      // 1) constructors, destructors, operators, etc.
      class CommandLineArgument : public Debug{
        
        public:
          
          CommandLineArgument(){
            boolValue_ = NULL;
            intValue_ = NULL;
            intValueList_ = NULL;
            doubleValue_ = NULL;
            doubleValueList_ = NULL;
            stringValue_ = NULL;
            stringValueList_ = NULL;
            isSet_ = false;
          };
                  
          int print(stringstream &s) const{
            s << "[CommandLine]   ";
            if((isAnOption_)||(isOptional_)){
              s << "[";
            }  
            
            s << "-" << key_;
            
            if(isAnOption_)
              s << ":";
            
            s << " ";
            if(!isAnOption_){
              s << "<";
              if((stringValueList_)||(intValueList_)||(doubleValueList_)){
                s << "{";
              }
            }
            
            if(description_.length()){
              s << description_;
            }
            else{
              s << "no description";
            }
       
            if((stringValue_)||(intValue_)||(doubleValue_)){
              s << " (default: ";
              if(stringValue_)
                s << "`" << *stringValue_ << "'";
              if(intValue_)
                s << *intValue_;
              if(doubleValue_)
                s << *doubleValue_;
              s << ")";
            }
            
            if(isAnOption_){
              s << " (default: " << *boolValue_ << ")";
            }

            if(!isAnOption_){
              if((stringValueList_)||(intValueList_)||(doubleValueList_)){
                s << "}";
              }
              s << ">";
            }
            
            if((isAnOption_)||(isOptional_)){
              s << "]";
            }
            s << endl;
            
            return 0;
          };
         
          bool              isOptional_, isAnOption_, isSet_;
          bool              *boolValue_;
          int               *intValue_;
          double            *doubleValue_;
          string            *stringValue_;
          vector<int>       *intValueList_;
          vector<double>    *doubleValueList_;
          vector<string>    *stringValueList_;
          
          string            key_;
          string            description_;
      };
      
      CommandLineParser(){
        setArgument("d", &(ttk::globalDebugLevel_),
          "Global debug level", true);
        setArgument("t", &threadNumber_, "Thread number", true);
      };
      
      ~CommandLineParser(){};
      
      // 2) functions
      int parse(int argc, char **argv){
       
        for(int i = 0; i < argc; i++){
          
          if((string(argv[i]) == "-h")||(string(argv[i]) == "--help")){
            printUsage(argv[0]);
          }
          
          for(int j = 0; j < (int) arguments_.size(); j++){
            
            if(!arguments_[j].isAnOption_){
              if((string(argv[i]) == "-" + arguments_[j].key_)
                &&(i + 1 < argc)){
              
                if(arguments_[j].stringValue_){
                  // let's process a string
                  (*arguments_[j].stringValue_) = argv[i + 1];
                  arguments_[j].isSet_ = true;
                }
                else if(arguments_[j].stringValueList_){
                  arguments_[j].stringValueList_->push_back(argv[i + 1]);
                  arguments_[j].isSet_ = true;
                }
                else if(arguments_[j].intValue_){
                  stringstream s(argv[i + 1]);
                  s >> *(arguments_[j].intValue_);
                  arguments_[j].isSet_ = true;
                }
                else if(arguments_[j].intValueList_){
                  stringstream s(argv[i + 1]);
                  arguments_[j].intValueList_->resize(
                    arguments_[j].intValueList_->size() + 1);
                  s >> arguments_[j].intValueList_->back();
                  arguments_[j].isSet_ = true;
                }
                else if(arguments_[j].doubleValue_){
                  stringstream s(argv[i + 1]);
                  s >> *(arguments_[j].doubleValue_);
                  arguments_[j].isSet_ = true;
                }
                else if(arguments_[j].doubleValueList_){
                  stringstream s(argv[i + 1]);
                  arguments_[j].doubleValueList_->resize(
                    arguments_[j].doubleValueList_->size() + 1);
                  s >> arguments_[j].doubleValueList_->back();
                  arguments_[j].isSet_ = true;
                }
              }
            }
            else{
              if(string(argv[i]) == "-" + arguments_[j].key_){
                *(arguments_[j].boolValue_) = !(*(arguments_[j].boolValue_));
              }
            }
          }
        }
        
        // check all the necessary arguments have been provided
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(!arguments_[i].isOptional_){
            if(!arguments_[i].isSet_){
              stringstream msg;
              msg << "[CommandLine] Missing mandatory argument:" << endl;
              arguments_[i].print(msg);
              dMsg(cerr, msg.str(), 1);
              printUsage(argv[0]);
            }
          }
        }
        
        return 0;
      };
      
      int printArgs(ostream &o = cout) const{
        
        o << "[CommandLine] Options and arguments:" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          o << "[CommandLine]   -" << arguments_[i].key_;
          o << ": ";
          
          if(arguments_[i].isAnOption_){
            if(arguments_[i].boolValue_){
              if(*(arguments_[i].boolValue_))
                o << "true";
              else
                o << "false";
            }
            else{
              o << "(not set)";
            }
          }
          else if(arguments_[i].stringValue_){
            if(arguments_[i].isSet_){
              o << *(arguments_[i].stringValue_);
            }
            else{
              o << "(not set)";
            }
          }
          else if(arguments_[i].stringValueList_){
            if(!arguments_[i].isSet_){
              o << "(not set)";
            }
            else{
              for(int j = 0; 
                j < (int) arguments_[i].stringValueList_->size(); j++){
                o << (*(arguments_[i].stringValueList_))[j] << " ";
              }
            }
          }
          else if(arguments_[i].intValue_){
            if(!arguments_[i].isSet_){
              o << "(not set)";
            }
            else{
              o << *(arguments_[i].intValue_);
            }
          }
          else if(arguments_[i].intValueList_){
            if(!arguments_[i].isSet_){
              o << "(not set)";
            }
            else{
              for(int j = 0; 
                j < (int) arguments_[i].intValueList_->size(); j++){
                o << (*(arguments_[i].intValueList_))[j] << " ";
              }
            }
          }
          else if(arguments_[i].doubleValue_){
            if(!arguments_[i].isSet_){
              o << "(not set)";
            }
            else{
              o << *(arguments_[i].doubleValue_);
            }
          }
          else if(arguments_[i].doubleValueList_){
            if(!arguments_[i].isSet_){
              o << "(not set)";
            }
            else{
              for(int j = 0; 
                j < (int) arguments_[i].doubleValueList_->size(); j++){
                o << (*(arguments_[i].doubleValueList_))[j] << " ";
              }
            }
          } 
          
          o << endl;
        }
        
        return 0;
      }
      
      int printUsage(const string &binPath) const {
        
        stringstream msg;
        msg << "[CommandLine]" << endl;
        msg << "[CommandLine] Usage:" << endl;
        msg << "[CommandLine]   " << binPath << endl;
        msg << "[CommandLine] Argument(s):" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(!arguments_[i].isAnOption_){
            arguments_[i].print(msg);
          }
        }
        msg << "[CommandLine] Option(s):" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(arguments_[i].isAnOption_){
            arguments_[i].print(msg);
          }
        }
        
        dMsg(cerr, msg.str(), 1);
        
        exit(0);
        return 0;
      };
      
      int setOption(const string &key, bool *value,
        const string &description = ""){
        
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
      
      int setArgument(const string &key, double *value, 
        const string &description = "",
        const bool &optional = false){
        
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
      
      int setArgument(const string &key, vector<double> *value, 
        const string &description = "",
        const bool &optional = false){
     
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
      
      inline int setArgument(const string &key, int *value, 
        const string &description = "",
        const bool &optional = false){
        
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
      
      int setArgument(const string &key, vector<int> *value, 
        const string &description = "",
        const bool &optional = false){
     
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
      
      int setArgument(const string &key, string *value, 
        const string &description = "",
        const bool &optional = false){
        
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
      
      int setArgument(const string &key, vector<string> *value, 
        const string &description = "",
        const bool &optional = false){
     
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
      
      vector<CommandLineArgument>
                        arguments_;
  };
}

#endif
