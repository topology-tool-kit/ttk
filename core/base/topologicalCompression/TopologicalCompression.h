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
#include <TopologicalSimplification.h>
#include <Triangulation.h>

// std
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stack>
#include <type_traits>

namespace ttk {

  enum class CompressionType { PersistenceDiagram = 0, Other = 1 };

  class TopologicalCompression : virtual public Debug {

  public:
    // Base code methods.
    TopologicalCompression();

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p inputOffsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    template <class dataType,
              typename triangulationType = AbstractTriangulation>
    int execute(const dataType *const inputData,
                const SimplexId *const inputOffsets,
                dataType *outputData,
                const triangulationType &triangulation);

    // Persistence compression methods.
    template <class dataType, typename triangulationType>
    int computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &JTPairs,
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &STPairs,
      const dataType *const inputScalars_,
      const SimplexId *const inputOffsets,
      const triangulationType &triangulation);
    template <typename dataType, typename triangulationType>
    int compressForPersistenceDiagram(int vertexNumber,
                                      const dataType *const inputData,
                                      const SimplexId *const inputOffset,
                                      dataType *outputData,
                                      const double &tol,
                                      const triangulationType &triangulation);

    // Other compression methods.
    template <typename dataType>
    int computeOther();

    template <typename dataType>
    int compressForOther(int vertexNumber,
                         const dataType *const inputData,
                         const SimplexId *const inputOffsets,
                         dataType *outputData,
                         const double &tol);

    // Getters and setters.
    inline void setCompressionType(int compressionType) {
      compressionType_ = compressionType;
    }
    inline void setSQ(std::string sqMethod) {
      SQMethod = sqMethod;
    }
    inline void setZFPOnly(bool z) {
      ZFPOnly = z;
    }
    inline void setSubdivide(bool b) {
      Subdivide = b;
    }
    inline void setMaximumError(double maximumError) {
      MaximumError = maximumError;
    }
    inline void setTolerance(const double data) {
      Tolerance = data;
    }
    inline void
      setUseTopologicalSimplification(bool useTopologicalSimplification) {
      UseTopologicalSimplification = useTopologicalSimplification;
    }
    inline void setFileName(char *fn) {
      fileName = fn;
    }
    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionVertexNeighbors();
        topologicalSimplification.preconditionTriangulation(triangulation);
        ftmTreePP.preconditionTriangulation(triangulation, false);
      }
    }

    inline int getNbVertices() {
      return NbVertices;
    }
    inline int getNbSegments() {
      return NbSegments;
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
      return SQMethodInt;
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
      return Tolerance;
    }
    inline double getZFPBitBudget() {
      return ZFPBitBudget;
    }
    inline bool getZFPOnly() {
      return ZFPOnly;
    }

    inline const std::vector<char> &getDataArrayName() const {
      return dataArrayName_;
    }
    inline std::vector<double> &getDecompressedData() {
      return decompressedData_;
    }
    inline std::vector<SimplexId> &getDecompressedOffsets() {
      return decompressedOffsets_;
    }
    inline std::vector<SimplexId> &getCompressedOffsets() {
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

    int ReadCompactSegmentation(FILE *fm,
                                std::vector<int> &segmentation,
                                int &numberOfVertices,
                                int &numberOfSegments);
    int ReadPersistenceIndex(
      FILE *fm,
      std::vector<std::tuple<double, int>> &mappings,
      std::vector<std::tuple<double, int>> &mappingsSortedPerValue,
      std::vector<std::tuple<int, double, int>> &constraints,
      double &min,
      double &max,
      int &nbConstraints);
    template <typename dataType>
    int ReadMetaData(FILE *fm);
    template <typename dataType, typename triangulationType>
    int ReadFromFile(FILE *fm, const triangulationType &triangulation);

    int WriteCompactSegmentation(FILE *fm,
                                 const std::vector<int> &segmentation,
                                 int numberOfVertices,
                                 int numberOfSegments);
    int WritePersistenceIndex(
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
    void CropIntervals(
      std::vector<std::tuple<dataType, int>> &mappings,
      std::vector<std::tuple<dataType, int>> &mappingsSortedPerValue,
      double min,
      double max,
      int vertexNumber,
      double *array,
      std::vector<int> &Seg);

    // API management.

#ifdef TTK_ENABLE_ZFP
    int CompressWithZFP(FILE *file,
                        bool decompress,
                        std::vector<double> &array,
                        int nx,
                        int ny,
                        int nz,
                        double rate);
#endif

#ifdef TTK_ENABLE_ZLIB
    unsigned long GetZlibDestLen(const unsigned long sourceLen);
    void CompressWithZlib(bool decompress,
                          unsigned char *dest,
                          unsigned long &destLen,
                          const unsigned char *const source,
                          const unsigned long sourceLen);
#endif

  private:
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
    template <typename dataType, typename triangulationType>
    int ReadPersistenceGeometry(FILE *fm,
                                const triangulationType &triangulation);
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

    template <typename dataType, typename triangulationType>
    int PerformSimplification(
      const std::vector<std::tuple<int, double, int>> &constraints,
      int nbConstraints,
      int vertexNumber,
      double *array,
      const triangulationType &triangulation);

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
    TopologicalSimplification topologicalSimplification{};
    ftm::FTMTreePP ftmTreePP;

    // Parameters
    int compressionType_{};
    bool ZFPOnly{false};
    double ZFPBitBudget{};
    int SQMethodInt{};

    double Tolerance{10};
    double MaximumError{10};
    int CompressionType{static_cast<int>(CompressionType::PersistenceDiagram)};
    std::string SQMethod{};
    bool Subdivide{false};
    bool UseTopologicalSimplification{true};

    int dataScalarType_{};
    int dataExtent_[6];
    double dataSpacing_[3];
    double dataOrigin_[3];

    std::vector<char> dataArrayName_{};

    // Persistence compression.
    std::vector<int> segmentation_{};
    std::vector<std::tuple<double, int>> mapping_{};
    std::vector<std::tuple<int, double, int>> criticalConstraints_{};

    // IO.
    int NbVertices{0};
    int NbSegments{0};
    int rawFileLength{0};
    std::vector<double> decompressedData_{};
    std::vector<SimplexId> decompressedOffsets_{};
    std::vector<SimplexId> compressedOffsets_{};
    int vertexNumberRead_{};
    char *fileName{};

    // Char array that identifies the file format.
    const char *magicBytes_{"TTKCompressedFileFormat"};
    // Current version of the file format. To be incremented at every
    // breaking change to keep backward compatibility.
    const unsigned long formatVersion_{1};
  };

} // namespace ttk

#include <OtherCompression.h>
#include <PersistenceDiagramCompression.h>

template <class dataType, typename triangulationType>
int ttk::TopologicalCompression::execute(
  const dataType *const inputData,
  const SimplexId *const inputOffsets,
  dataType *outputData,
  const triangulationType &triangulation) {
  this->printMsg("Starting compression...");

// check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(inputData == nullptr)
    return -2;
  if(outputData == nullptr)
    return -3;
// if (tol < 0 || tol > 100) return -4;
#endif

  int vertexNumber = triangulation.getNumberOfVertices();

  int res = 0;
  if(compressionType_ == (int)ttk::CompressionType::PersistenceDiagram)
    compressForPersistenceDiagram(vertexNumber, inputData, inputOffsets,
                                  outputData, Tolerance, triangulation);
  else if(compressionType_ == (int)ttk::CompressionType::Other)
    compressForOther(
      vertexNumber, inputData, inputOffsets, outputData, Tolerance);

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
  NbVertices = numberOfVertices;

  int totalSize = usePersistence
                    ? ComputeTotalSizeForPersistenceDiagram<double>(
                      getMapping(), getCriticalConstraints(), zfpOnly,
                      getNbSegments(), getNbVertices(), zfpBitBudget)
                  : useOther ? ComputeTotalSizeForOther<double>()
                             : 0;

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
    this->printMsg("Geometry successfully written to buffer.");
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
  auto sourceLen = static_cast<unsigned long>(rawFileLength);
  const auto source = reinterpret_cast<unsigned char *>(buf);
  auto destLen = GetZlibDestLen(sourceLen);
  std::vector<unsigned char> ddest(destLen);
  CompressWithZlib(false, ddest.data(), destLen, source, sourceLen);
  this->printMsg("Data successfully compressed.");

  // [fm->fp] Copy fm to fp.
  Write(fp, destLen); // Compressed size...
  Write(fp, sourceLen);
  WriteByteArray(fp, ddest.data(), destLen);
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
  int sqType = (strcmp(sq, "") == 0)                            ? 0
               : (strcmp(sq, "r") == 0 || strcmp(sq, "R") == 0) ? 1
               : (strcmp(sq, "d") == 0 || strcmp(sq, "D") == 0) ? 2
                                                                : 3;

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
  // (explicit call to unsigned long variant for MSVC compatibility)
  Write<unsigned long>(fp, dataArrayName.size());

  // 7. Array name (as unsigned chars)
  WriteByteArray(fp, dataArrayName.c_str(), dataArrayName.size());

  this->printMsg("Metadata successfully written.");

  return 0;
}

template <typename T, typename triangulationType>
int ttk::TopologicalCompression::ReadFromFile(
  FILE *fp, const triangulationType &triangulation) {
  // [fp->] Read headers.

  this->printMsg("Successfully read metadata.");

  if(ZFPOnly && (ZFPBitBudget > 64 || ZFPBitBudget < 1)) {
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
    auto sl = Read<unsigned long>(fp); // Compressed size...
    auto dl = Read<unsigned long>(fp); // Uncompressed size...

    destLen = dl;
    std::vector<unsigned char> ssource(sl);
    ReadByteArray(fp, ssource.data(), sl);
    this->printMsg("Successfully read compressed data.");

    // [ff->fm] Decompress data.
    ddest.resize(destLen);
    dest = ddest.data();
    CompressWithZlib(true, dest, destLen, ssource.data(), sl);
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
  if(!(ZFPOnly)) {
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
    status = ReadPersistenceGeometry<double>(fm, triangulation);
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
  ZFPOnly = Read<bool>(fm);

  // 0. SQ type
  SQMethodInt = Read<int>(fm);

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
  Tolerance = Read<double>(fm);

  // 5. Lossy compressor ratio
  ZFPBitBudget = Read<double>(fm);

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
