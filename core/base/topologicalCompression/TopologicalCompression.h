/// \ingroup baseCode
/// \class ttk::TopologicalCompression
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 20/04/2017
///
/// \brief TTK %topologicalCompression processing package.
///
/// %TopologicalCompression is a TTK processing package that takes a scalar
/// field on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkTopologicalCompression.cpp %for a usage example.

#pragma once

// base code includes
#include <FTMTreePP.h>
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>
#include <Triangulation.h>

// std

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stack>
#include <type_traits>

#ifdef TTK_ENABLE_ZLIB
#include <zlib.h>
#endif

#ifdef TTK_ENABLE_ZFP
#ifndef __cplusplus
#define __cplusplus 201112L
#endif
#include <climits>
#include <cstdio>
#include <zfp.h>
#include <zfp/macros.h>
#endif

namespace ttk {

  enum class CompressionType { PersistenceDiagram = 0, Other = 1 };

  class TopologicalCompression : virtual public Debug {

  public:
    // Base code methods.
    TopologicalCompression();

    template <class dataType>
    int execute(const double &tolerance);

    // Persistence compression methods.
    template <class dataType>
    int computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &JTPairs,
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &STPairs,
      dataType *inputScalars_,
      SimplexId *inputOffsets);
    template <typename dataType>
    int compressForPersistenceDiagram(int vertexNumber,
                                      dataType *inputData,
                                      dataType *outputData,
                                      const double &tol);

    // Other compression methods.
    template <typename dataType>
    int computeOther();

    template <typename dataType>
    int compressForOther(int vertexNumber,
                         dataType *inputData,
                         dataType *outputData,
                         const double &tol);

    // Getters and setters.
    inline void setInputDataPointer(void *data) {
      inputData_ = data;
    }
    inline void setOutputDataPointer(void *data) {
      outputData_ = data;
    }
    inline void setCompressionType(int compressionType) {
      compressionType_ = compressionType;
    }
    inline void setSQ(std::string sqMethod) {
      sqMethod_ = sqMethod;
    }
    inline void setZFPOnly(bool z) {
      zfpOnly_ = z;
    }
    inline void setSubdivide(bool dontSubdivide) {
      dontSubdivide_ = dontSubdivide;
    }
    inline void setMaximumError(double maximumError) {
      maximumError_ = maximumError;
    }
    inline void
      setUseTopologicalSimplification(bool useTopologicalSimplification) {
      useTopologicalSimplification_ = useTopologicalSimplification;
    }
    inline void setFileName(char *fn) {
      fileName = fn;
    }
    inline void setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_)
        triangulation_->preconditionVertexNeighbors();
    }

    inline int getNbVertices() {
      return nbVertices;
    }
    inline int getNbSegments() {
      return nbSegments;
    }
    inline std::vector<int> &getSegmentation() {
      return segmentation_;
    }
    inline std::vector<std::tuple<double, int>> &getMapping() {
      return mapping_;
    }
    inline std::vector<std::tuple<int, double, int>> &getCriticalConstraints() {
      return criticalConstraints_;
    }

    inline int getCompressionType() {
      return compressionType_;
    }
    inline int getSQMethod() {
      return sqMethodInt_;
    }
    inline int getDataScalarType() {
      return dataScalarType_;
    }
    inline int *getDataExtent() {
      return dataExtent_;
    }
    inline double *getDataSpacing() {
      return dataSpacing_;
    }
    inline double *getDataOrigin() {
      return dataOrigin_;
    }
    inline double getTolerance() {
      return tolerance_;
    }
    inline double getZFPBitBudget() {
      return zfpBitBudget_;
    }
    inline bool getZFPOnly() {
      return zfpOnly_;
    }

    inline const std::vector<char> &getDataArrayName() const {
      return dataArrayName_;
    }
    inline std::vector<double> &getDecompressedData() {
      return decompressedData_;
    }
    inline std::vector<int> &getDecompressedOffsets() {
      return decompressedOffsets_;
    }
    inline std::vector<int> &getCompressedOffsets() {
      return compressedOffsets_;
    }

    // IO management.
    static unsigned int log2(int val);
    inline static bool cmp(const std::tuple<double, int> &a,
                           const std::tuple<double, int> &b) {
      return std::get<1>(a) > std::get<1>(b);
    };
    inline static bool cmp2(const std::tuple<double, int> &a,
                            const std::tuple<double, int> &b) {
      return std::get<1>(a) < std::get<1>(b);
    };

    template <typename T>
    static T Read(FILE *fm) {
      static_assert(std::is_same<T, bool>() || std::is_same<T, int>()
                      || std::is_same<T, unsigned long>()
                      || std::is_same<T, double>(),
                    "Function should have only those types");
      T ret;
      const auto status = std::fread(&ret, sizeof(T), 1, fm);
      if(status == 0) {
        ttk::Debug dbg;
        dbg.printErr("Error reading " + std::string(typeid(T).name()) + "!");
      }
      return ret;
    }
    template <typename T>
    static void ReadByteArray(FILE *fm, T *buffer, size_t length) {
      const auto status = std::fread(buffer, sizeof(T), length, fm);
      if(status == 0) {
        ttk::Debug dbg;
        dbg.printErr("Error reading " + std::string(typeid(T).name())
                     + "array!");
      }
    }
    template <typename T>
    static void Write(FILE *fm, T data) {
      static_assert(std::is_same<T, bool>() || std::is_same<T, int>()
                      || std::is_same<T, unsigned long>()
                      || std::is_same<T, double>(),
                    "Function should have only those types");
      auto status = std::fwrite(&data, sizeof(T), 1, fm);
      if(status == 0) {
        ttk::Debug dbg;
        dbg.printErr("Error writing " + std::string(typeid(T).name()) + "!");
      }
    }
    template <typename T>
    static void WriteByteArray(FILE *fm, const T *buffer, size_t length) {
      const auto status = std::fwrite(buffer, sizeof(T), length, fm);
      if(status == 0) {
        ttk::Debug dbg;
        dbg.printErr("Error writing " + std::string(typeid(T).name())
                     + "array!");
      }
    }

    static int ReadCompactSegmentation(FILE *fm,
                                       std::vector<int> &segmentation,
                                       int &numberOfVertices,
                                       int &numberOfSegments);
    static int ReadPersistenceIndex(
      FILE *fm,
      std::vector<std::tuple<double, int>> &mappings,
      std::vector<std::tuple<double, int>> &mappingsSortedPerValue,
      std::vector<std::tuple<int, double, int>> &constraints,
      double &min,
      double &max,
      int &nbConstraints);
    template <typename dataType>
    int ReadMetaData(FILE *fm);
    template <typename dataType>
    int ReadFromFile(FILE *fm);

    static int WriteCompactSegmentation(FILE *fm,
                                        const std::vector<int> &segmentation,
                                        int numberOfVertices,
                                        int numberOfSegments);
    static int WritePersistenceIndex(
      FILE *fm,
      std::vector<std::tuple<double, int>> &mapping,
      std::vector<std::tuple<int, double, int>> &constraints);
    template <typename T>
    int WriteMetaData(FILE *fp,
                      int compressionType,
                      bool zfpOnly,
                      const char *sqMethod,
                      int dataType,
                      int *dataExtent,
                      double *dataSpacing,
                      double *dataOrigin,
                      double tolerance,
                      double zfpBitBudget,
                      const std::string &dataArrayName);
    template <typename T>
    int WriteToFile(FILE *fp,
                    int compressionType,
                    bool zfpOnly,
                    const char *sqMethod,
                    int dataType,
                    int *dataExtent,
                    double *dataSpacing,
                    double *dataOrigin,
                    double *data,
                    double tolerance,
                    double zfpBitBudget,
                    const std::string &dataArrayName);

    template <typename dataType>
    static void CropIntervals(
      std::vector<std::tuple<dataType, int>> &mappings,
      std::vector<std::tuple<dataType, int>> &mappingsSortedPerValue,
      double min,
      double max,
      int vertexNumber,
      double *array,
      std::vector<int> &Seg);

    // API management.

#ifdef TTK_ENABLE_ZFP
    static int CompressWithZFP(FILE *file,
                               bool decompress,
                               std::vector<double> &array,
                               int nx,
                               int ny,
                               int nz,
                               double rate);
#endif

#ifdef TTK_ENABLE_ZLIB
    static void CompressWithZlib(bool decompress,
                                 Bytef *dest,
                                 uLongf *destLen,
                                 const Bytef *source,
                                 uLong sourceLen);
#endif

  private:
#ifdef TTK_ENABLE_ZFP
    static int compressZFPInternal(double *array,
                                   int nx,
                                   int ny,
                                   int nz,
                                   double rate,
                                   bool decompress,
                                   FILE *file);
#endif

    // Internal read/write.

    template <typename dataType>
    int ComputeTotalSizeForOther();
    template <typename dataType>
    int ComputeTotalSizeForPersistenceDiagram(
      std::vector<std::tuple<double, int>> &mapping,
      std::vector<std::tuple<int, double, int>> &criticalConstraints,
      bool zfpOnly,
      int nbSegments,
      int nbVertices,
      double zfpBitBudget);

    template <typename dataType>
    int ReadPersistenceTopology(FILE *fm);
    template <typename dataType>
    int ReadOtherTopology(FILE *fm);
    template <typename dataType>
    int ReadPersistenceGeometry(FILE *fm);
    template <typename dataType>
    int ReadOtherGeometry(FILE *fm);

    template <typename dataType>
    int WritePersistenceTopology(FILE *fm);
    template <typename dataType>
    int WriteOtherTopology(FILE *fm);
    template <typename dataType>
    int WritePersistenceGeometry(FILE *fm,
                                 int *dataExtent,
                                 bool zfpOnly,
                                 double zfpBitBudget,
                                 double *toCompress);
    template <typename dataType>
    int WriteOtherGeometry(FILE *fm);

    template <typename dataType>
    int PerformSimplification(
      const std::vector<std::tuple<int, double, int>> &constraints,
      int nbConstraints,
      int vertexNumber,
      double *array);

    // Numeric management.

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    template <typename type>
    static type abs_diff(const type var1, const type var2) {
      return (var1 > var2) ? var1 - var2 : var2 - var1;
    }

  protected:
    // General.
    void *inputData_{};
    void *outputData_{};
    Triangulation *triangulation_{};
    TopologicalSimplification topologicalSimplification{};

    // Parameters
    int compressionType_{};
    bool zfpOnly_{};
    int sqMethodInt_{};
    std::string sqMethod_{""};
    int dataScalarType_{};
    int dataExtent_[6];
    double dataSpacing_[3];
    double dataOrigin_[3];
    double tolerance_{};
    double maximumError_{};
    double zfpBitBudget_{};
    bool dontSubdivide_{};
    bool useTopologicalSimplification_{};
    std::vector<char> dataArrayName_{};

    // Persistence compression.
    std::vector<int> segmentation_{};
    std::vector<std::tuple<double, int>> mapping_{};
    std::vector<std::tuple<int, double, int>> criticalConstraints_{};

    // IO.
    int nbVertices{0};
    int nbSegments{0};
    int rawFileLength{0};
    std::vector<double> decompressedData_{};
    std::vector<int> decompressedOffsets_{};
    std::vector<int> compressedOffsets_{};
    int vertexNumberRead_{};
    char *fileName{};

    // Char array that identifies the file format.
    static const char *magicBytes_;
    // Current version of the file format. To be incremented at every
    // breaking change to keep backward compatibility.
    static const unsigned long formatVersion_;
  };

} // namespace ttk

#include <OtherCompression.h>
#include <PersistenceDiagramCompression.h>

template <class dataType>
int ttk::TopologicalCompression::execute(const double &tol) {
  this->printMsg("Starting compression...");

// check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
  if(!outputData_)
    return -3;
// if (tol < 0 || tol > 100) return -4;
#endif

  auto *outputData = (dataType *)outputData_;
  auto *inputData = (dataType *)inputData_;
  int vertexNumber = triangulation_->getNumberOfVertices();

  int res = 0;
  if(compressionType_ == (int)ttk::CompressionType::PersistenceDiagram)
    compressForPersistenceDiagram<dataType>(
      vertexNumber, inputData, outputData, tol);
  else if(compressionType_ == (int)ttk::CompressionType::Other)
    compressForOther<dataType>(vertexNumber, inputData, outputData, tol);

  return res;
}

template <typename T>
int ttk::TopologicalCompression::WriteToFile(FILE *fp,
                                             int compressionType,
                                             bool zfpOnly,
                                             const char *sqMethod,
                                             int dataType,
                                             int *dataExtent,
                                             double *dataSpacing,
                                             double *dataOrigin,
                                             double *data,
                                             double tolerance,
                                             double zfpBitBudget,
                                             const std::string &dataArrayName) {
  // [->fp] Write metadata.
  WriteMetaData<double>(fp, compressionType, zfpOnly, sqMethod, dataType,
                        dataExtent, dataSpacing, dataOrigin, tolerance,
                        zfpBitBudget, dataArrayName);

#ifdef TTK_ENABLE_ZLIB
  Write(fp, true);
#else
  Write(fp, false);
#endif

  bool usePersistence
    = compressionType == (int)ttk::CompressionType::PersistenceDiagram;
  bool useOther = compressionType == (int)ttk::CompressionType::Other;

  int numberOfVertices = 1;
  for(int i = 0; i < 3; ++i)
    numberOfVertices *= (1 + dataExtent[2 * i + 1] - dataExtent[2 * i]);
  nbVertices = numberOfVertices;

  int totalSize = usePersistence
                    ? ComputeTotalSizeForPersistenceDiagram<double>(
                      getMapping(), getCriticalConstraints(), zfpOnly,
                      getNbSegments(), getNbVertices(), zfpBitBudget)
                    : useOther ? ComputeTotalSizeForOther<double>() : 0;

  std::vector<char> bbuf(totalSize);
  char *buf = bbuf.data();
  size_t len = (size_t)totalSize;

  // #ifndef _MSC_VER
  // FILE *fm = fmemopen(buf, len, "r+");
  // #else
  const std::string s = fileName + std::string(".temp");
  const char *ffn = s.c_str();
  FILE *fm = fopen(ffn, "wb");
  // #endif

  // [->fm] Encode, lossless compress and write topology.
  if(!(zfpOnly)) {
    if(usePersistence)
      WritePersistenceTopology<double>(fm);
    else if(useOther)
      WriteOtherTopology<double>(fm);
  }

  this->printMsg("Topology successfully written to buffer.");

  int status = 0;
  // [->fm] Write altered geometry.
  if(usePersistence)
    status = WritePersistenceGeometry<double>(
      fm, dataExtent, zfpOnly, zfpBitBudget, data);
  else if(useOther)
    status = WriteOtherGeometry<double>(fm);

  fclose(fm); // !Close stream to write changes!
  // #ifdef _MSC_VER
  fm = fopen(ffn, "rb");
  int ret = fread(buf, len, sizeof(char), fm);
  fclose(fm);
  remove(ffn);
  // #endif

  if(status == 0) {
    this->printMsg(" Geometry successfully written to buffer.");
  } else {
    this->printErr("Geometry was not successfully written to buffer.");
    fflush(fp);
    fclose(fp);
    return -1;
  }

  // Check computed size vs read size.
  if(totalSize < rawFileLength) {
    this->printErr("Invalid total size (" + std::to_string(totalSize) + " vs "
                   + std::to_string(rawFileLength) + ").");
  }

#ifdef TTK_ENABLE_ZLIB
  // [fm->ff] Compress fm.
  uLong sourceLen = (uLong)rawFileLength;
  Bytef *source = reinterpret_cast<unsigned char *>(buf);
  uLongf destLen = compressBound(sourceLen);
  std::vector<Bytef> ddest(destLen);
  Bytef *dest = ddest.data();
  CompressWithZlib(false, dest, &destLen, source, sourceLen);
  this->printMsg("Data successfully compressed.");

  // [fm->fp] Copy fm to fp.
  Write(fp, destLen); // Compressed size...
  Write(fp, sourceLen);
  WriteByteArray(fp, dest, destLen);
  this->printMsg("Data successfully written to filesystem.");

#else
  this->printMsg("ZLIB not found, writing raw file.");
  unsigned char *source = reinterpret_cast<unsigned char *>(buf);
  unsigned long sourceLen = (unsigned long)rawFileLength;
  unsigned long destLen = (unsigned long)rawFileLength;

  Write(fp, destLen); // Compressed size...
  Write(fp, sourceLen);
  WriteByteArray(fp, source, destLen);
#endif

  fflush(fp);
  fclose(fp);

  return ret;
}

template <typename T>
int ttk::TopologicalCompression::WriteMetaData(
  FILE *fp,
  int compressionType,
  bool zfpOnly,
  const char *sqMethod,
  int dataType,
  int *dataExtent,
  double *dataSpacing,
  double *dataOrigin,
  double tolerance,
  double zfpBitBudget,
  const std::string &dataArrayName) {

  // -4. Magic bytes
  WriteByteArray(fp, magicBytes_, std::strlen(magicBytes_));

  // -3. File format version
  Write(fp, formatVersion_);

  // -2. Persistence, or Other
  Write(fp, compressionType);

  // -1. zfpOnly
  Write(fp, zfpOnly);

  // 0. SQ type
  const char *sq = sqMethod;
  int sqType = (strcmp(sq, "") == 0)
                 ? 0
                 : (strcmp(sq, "r") == 0 || strcmp(sq, "R") == 0)
                     ? 1
                     : (strcmp(sq, "d") == 0 || strcmp(sq, "D") == 0) ? 2 : 3;

  Write(fp, sqType);

  // 1. DataType
  Write(fp, dataType);

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    Write(fp, dataExtent[i]);

  for(int i = 0; i < 3; ++i)
    Write(fp, dataSpacing[i]);

  for(int i = 0; i < 3; ++i)
    Write(fp, dataOrigin[i]);

  // 4. Tolerance
  Write(fp, tolerance);

  // 5. ZFP ratio
  Write(fp, zfpBitBudget);

  // 6. Length of array name
  Write(fp, dataArrayName.size());

  // 7. Array name (as unsigned chars)
  WriteByteArray(fp, dataArrayName.c_str(), dataArrayName.size());

  this->printMsg("Metadata successfully written.");

  return 0;
}

template <typename T>
int ttk::TopologicalCompression::ReadFromFile(FILE *fp) {
  // [fp->] Read headers.

  this->printMsg("Successfully read metadata.");

  if(zfpOnly_ && (zfpBitBudget_ > 64 || zfpBitBudget_ < 1)) {
    this->printMsg("Wrong ZFP bit budget for ZFP-only use.");
    return -4;
  }

  bool useZlib = Read<bool>(fp);
  unsigned char *dest;
  std::vector<unsigned char> ddest;
  unsigned long destLen;

#ifdef TTK_ENABLE_ZLIB
  if(useZlib) {
    // [fp->ff] Read compressed data.
    uLongf sl = Read<unsigned long>(fp); // Compressed size...
    uLongf dl = Read<unsigned long>(fp); // Uncompressed size...

    unsigned long sourceLen = (uLongf)sl;
    destLen = dl;
    std::vector<Bytef> ssource(sl);
    Bytef *source = ssource.data();
    ReadByteArray(fp, source, sl);
    this->printMsg("Successfully read compressed data.");

    // [ff->fm] Decompress data.
    ddest.resize(destLen);
    dest = ddest.data();
    CompressWithZlib(true, dest, &destLen, source, sourceLen);
    this->printMsg("Successfully uncompressed data.");

  } else {
    this->printMsg("File was not compressed with ZLIB.");

    Read<unsigned long>(fp); // Compressed size...
    unsigned long dl = Read<unsigned long>(fp); // Uncompressed size...

    destLen = dl;
    ddest.resize(destLen);
    dest = ddest.data();
    ReadByteArray(fp, dest, destLen);
  }
#else
  if(useZlib) {
    this->printMsg(" File compressed but ZLIB not installed! Aborting.");
    return -4;

  } else {
    this->printMsg(" ZLIB not installed, but file was not compressed anyways.");

    Read<unsigned long>(fp); // Compressed size...
    unsigned long dl = Read<unsigned long>(fp); // Uncompressed size...

    destLen = dl;
    ddest.resize(destLen);
    dest = ddest.data();
    ReadByteArray(fp, dest, destLen);
  }
#endif

  // [fm->] Read data.
  char *buf = reinterpret_cast<char *>(dest);

  //#ifndef _MSC_VER
  // FILE *fm = fmemopen(buf, destLen, "r+");
  //#else
  const std::string s = fileName + std::string(".temp");
  const char *ffn = s.c_str();
  FILE *ftemp = fopen(ffn, "wb");
  fwrite(buf, destLen, sizeof(char), ftemp);
  fclose(ftemp);
  FILE *fm = fopen(ffn, "rb");
  //#endif

  // Do read topology.
  if(!(zfpOnly_)) {
    if(compressionType_ == (int)ttk::CompressionType::PersistenceDiagram)
      ReadPersistenceTopology<double>(fm);
    else if(compressionType_ == (int)ttk::CompressionType::Other)
      ReadOtherTopology<double>(fm);
  }

  this->printMsg("Successfully read topology.");

  // Get altered geometry.
  // Rebuild topologically consistent geometry.
  int status = 0;
  if(compressionType_ == (int)ttk::CompressionType::PersistenceDiagram)
    status = ReadPersistenceGeometry<double>(fm);
  else if(compressionType_ == (int)ttk::CompressionType::Other)
    status = ReadOtherGeometry<double>(fm);

  fclose(fm);
  // #ifdef _MSC_VER
  remove(ffn);
  // #endif
  fclose(fp);

  if(status == 0) {
    this->printMsg("Successfully read geometry.");
    this->printMsg("Successfully read file.");
  } else {
    this->printMsg("Failed to write (possibly ZFP)!");
    this->printMsg("File may be corrupted!");
  }

  return status;
}

template <typename T>
int ttk::TopologicalCompression::ReadMetaData(FILE *fm) {

  // -4. Magic bytes
  const auto magicBytesLen{std::strlen(magicBytes_)};
  std::vector<char> mBytes(magicBytesLen + 1);
  mBytes[magicBytesLen] = '\0'; // NULL-termination
  ReadByteArray(fm, mBytes.data(), magicBytesLen);

  // To deal with pre-v1 file format (without scalar field array name)
  bool hasMagicBytes = strcmp(mBytes.data(), magicBytes_) == 0;

  if(!hasMagicBytes) {
    this->printWrn("Could not find magic bytes in input file!");
    this->printWrn("File might be corrupted!");

    // rewind fm to beginning of file
    std::rewind(fm);
  }

  // -3. File format version
  unsigned long version = 0;
  if(hasMagicBytes) {
    version = Read<unsigned long>(fm);
  }

  // -2. Compression type.
  compressionType_ = Read<int>(fm);

  // -1. ZFP only type.
  zfpOnly_ = Read<bool>(fm);

  // 0. SQ type
  sqMethodInt_ = Read<int>(fm);

  // 1. DataType
  dataScalarType_ = Read<int>(fm);
  // DataScalarType = VTK_DOUBLE;

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    dataExtent_[i] = Read<int>(fm);

  for(int i = 0; i < 3; ++i)
    dataSpacing_[i] = Read<double>(fm);

  for(int i = 0; i < 3; ++i)
    dataOrigin_[i] = Read<double>(fm);

  // 4. Error tolerance (relative percentage)
  tolerance_ = Read<double>(fm);

  // 5. Lossy compressor ratio
  zfpBitBudget_ = Read<double>(fm);

  if(version == 0) {
    // Pre-v1 format has no scalar field array name
    return 0;
  }

  // 6. Length of array name
  size_t dataArrayNameLength = Read<unsigned long>(fm);

  // 7. Array name (as unsigned chars)
  dataArrayName_.resize(dataArrayNameLength + 1);
  dataArrayName_[dataArrayNameLength] = '\0'; // NULL-termination
  ReadByteArray(fm, dataArrayName_.data(), dataArrayNameLength);

  return 0;
}
