#include <Debug.h>
#include <Os.h>

#ifdef _WIN32

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <windows.h>

#include <ciso646>
#include <cwchar>
#include <direct.h>
#include <stdint.h>
#include <time.h>

#elif defined(__unix__) || defined(__APPLE__)

#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <iostream>
#include <sstream>

namespace ttk {

  int OsCall::getCurrentDirectory(std::string &directoryPath) {
#ifdef _WIN32
    directoryPath = _getcwd(NULL, 0);
#else
    directoryPath = getcwd(NULL, PATH_MAX);
#endif
    directoryPath += "/";

    return 0;
  }

  float OsCall::getMemoryInstantUsage() {
#ifdef __linux__
    // horrible hack since getrusage() doesn't seem to work well under
    // linux
    std::stringstream procFileName;
    procFileName << "/proc/" << getpid() << "/statm";

    std::ifstream procFile(procFileName.str().data(), std::ios::in);
    if(procFile) {
      float memoryUsage;
      procFile >> memoryUsage;
      procFile.close();
      return memoryUsage / 1024.0;
    }
#endif
    return 0;
  }

  int OsCall::getNumberOfCores() {
#ifdef TTK_ENABLE_OPENMP
    return omp_get_num_procs();
#endif
    return 1;
  }

  double OsCall::getTimeStamp() {
#ifdef _WIN32
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);

    LARGE_INTEGER temp;
    QueryPerformanceCounter(&temp);

    return (double)temp.QuadPart / frequency.QuadPart;
#endif

#ifdef __APPLE__
    struct timeval stamp;
    gettimeofday(&stamp, NULL);
    return (stamp.tv_sec * 1000000 + stamp.tv_usec) / 1000000.0;
#endif

#ifdef __unix__
    struct timeval stamp;
    gettimeofday(&stamp, NULL);
    return (stamp.tv_sec * 1000000 + stamp.tv_usec) / 1000000.0;
#endif
  }

  std::vector<std::string>
    OsCall::listFilesInDirectory(const std::string &directoryName,
                                 const std::string &extension) {

    std::vector<std::string> filesInDir;

#ifdef _WIN32

#ifdef UNICODE
    auto toWString = [](const std::string &str) {
      if(str.empty())
        return std::wstring();
      int wcharCount
        = MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, nullptr, 0);
      std::wstring wstr;
      wstr.resize(wcharCount);
      MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, &wstr[0], wcharCount);
      return wstr;
    };
    auto toString = [](const WCHAR wstr[]) {
      int charCount = WideCharToMultiByte(
        CP_UTF8, 0, wstr, -1, nullptr, 0, nullptr, nullptr);
      std::string str;
      str.resize(charCount);
      WideCharToMultiByte(
        CP_UTF8, 0, wstr, -1, &str[0], charCount, nullptr, nullptr);
      return str;
    };
#else
    auto toWString = [](const std::string &str) { return str; };
    auto toString = [](const char *c) { return std::string(c); };
#endif

    WIN32_FIND_DATA FindFileData;
    char *buffer;
    buffer = _getcwd(NULL, 0);
    if((buffer = _getcwd(NULL, 0)) == NULL)
      perror("_getcwd error");
    else {
      free(buffer);
    }

    HANDLE hFind
      = FindFirstFile(toWString(directoryName).c_str(), &FindFileData);
    if(hFind == INVALID_HANDLE_VALUE) {
      std::stringstream msg;
      msg << "[Os] Could not open directory `" << directoryName
          << "'. Error: " << GetLastError() << std::endl;
      Debug d;
      d.dMsg(std::cerr, msg.str(), 0);
    } else {
      const std::string filename = toString(FindFileData.cFileName);

      std::string entryExtension(filename);
      entryExtension
        = entryExtension.substr(entryExtension.find_last_of('.') + 1);

      if(entryExtension == extension)
        filesInDir.push_back(filename);
      std::string dir = directoryName;
      dir.resize(dir.size() - 1);
      while(FindNextFile(hFind, &FindFileData)) {
        if(extension.size()) {
          std::string entryExtension(filename);
          entryExtension
            = entryExtension.substr(entryExtension.find_last_of('.') + 1);
          if(entryExtension == extension)
            filesInDir.push_back(dir + filename);
        } else {
          if((filename != ".") && (filename != "..")) {
            filesInDir.push_back(directoryName + "/" + filename);
          }
        }
      }
    }
    FindClose(hFind);
#else
    DIR *d = opendir((directoryName + "/").data());
    if(!d) {
      std::stringstream msg;
      msg << "[Os] Could not open directory `" << directoryName << "'..."
          << std::endl;
      Debug dbg;
      dbg.dMsg(std::cerr, msg.str(), 0);
    } else {
      struct dirent *dirEntry;
      while((dirEntry = readdir(d)) != NULL) {
        if(extension.size()) {
          std::string entryExtension(dirEntry->d_name);
          entryExtension
            = entryExtension.substr(entryExtension.find_last_of('.') + 1);
          if(entryExtension == extension)
            filesInDir.push_back(directoryName + "/"
                                 + std::string(dirEntry->d_name));
        } else {
          if((std::string(dirEntry->d_name) != ".")
             && (std::string(dirEntry->d_name) != ".."))
            filesInDir.push_back(directoryName + "/"
                                 + std::string(dirEntry->d_name));
        }
      }
      closedir(d);
    }
#endif

    std::sort(filesInDir.begin(), filesInDir.end());

    return filesInDir;
  }

  int OsCall::mkDir(const std::string &directoryName) {

#ifdef _WIN32
    return _mkdir(directoryName.data());
#else
    return mkdir(directoryName.data(), 0777);
#endif
  }

  int OsCall::nearbyint(const double &x) {
    const double upperBound = ceil(x);
    const double lowerBound = floor(x);

    if(upperBound - x <= x - lowerBound)
      return (int)upperBound;
    else
      return (int)lowerBound;
  }

  int OsCall::rmDir(const std::string &directoryName) {

#ifdef _WIN32
    // NOTE:
    // the directory will be deleted with this call
    // only if it's empty...
    return _rmdir(directoryName.data());
#else
    std::stringstream cmd;
    cmd << "rm -R " << directoryName << " 2> /dev/null";
    return system(cmd.str().data());
#endif
  }

  int OsCall::rmFile(const std::string &fileName) {

    std::stringstream cmd;

#ifdef _WIN32
    cmd << "del";
#else
    cmd << "rm";
#endif

    cmd << " " << fileName;

#ifndef _WIN32
    cmd << " 2> /dev/null";
#endif

    return system(cmd.str().data());
  }

} // namespace ttk
