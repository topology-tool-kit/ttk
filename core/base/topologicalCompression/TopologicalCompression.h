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

#ifndef _TOPOLOGICALCOMPRESSION_H
#define _TOPOLOGICALCOMPRESSION_H

// base code includes

#include <DataTypes.h>
#include <Debug.h>
#include <FTMTreePP.h>
#include <PersistenceDiagram.h>
#include <TopologicalSimplification.h>
#include <Triangulation.h>
#include <Wrapper.h>

// std

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <string.h>

#ifdef TTK_ENABLE_ZLIB
#include <zlib.h>
#endif

#ifdef TTK_ENABLE_ZFP
#ifndef __cplusplus
#define __cplusplus 201112L
#endif
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zfp.h>
#include <zfp/macros.h>
#endif

namespace ttk {

  enum class CompressionType { PersistenceDiagram = 0, Other = 1 };

  class TopologicalCompression : public Debug {

  public:
    // Base code methods.
    TopologicalCompression();
    ~TopologicalCompression();
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
    inline int setInputDataPointer(void *data) {
      inputData_ = data;
      return 0;
    }

    inline int setOutputDataPointer(void *data) {
      outputData_ = data;
      return 0;
    }

    inline int setCompressionType(int compressionType) {
      compressionType_ = compressionType;
      return 0;
    }

    inline int setSQ(std::string sqMethod) {
      sqMethod_ = sqMethod;
      return 0;
    }

    inline int setZFPOnly(bool z) {
      zfpOnly_ = z;
      return 0;
    }

    inline int setSubdivide(bool dontSubdivide) {
      dontSubdivide_ = dontSubdivide;
      return 0;
    }

    inline int setMaximumError(double maximumError) {
      maximumError_ = maximumError;
      return 0;
    }

    inline int
      setUseTopologicalSimplification(bool useTopologicalSimplification) {
      useTopologicalSimplification_ = useTopologicalSimplification;
      return 0;
    }

    inline int setFileName(char *fn) {
      fileName = fn;
      return 0;
    }

    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_)
        triangulation_->preprocessVertexNeighbors();
      return 0;
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

    static bool ReadBool(FILE *fm);
    static int ReadInt(FILE *fm);
    static double ReadDouble(FILE *fm);
    static unsigned long ReadUnsignedLong(FILE *fm);
    static void
      ReadUnsignedCharArray(FILE *fm, unsigned char *buffer, size_t length);
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

    static void WriteBool(FILE *fm, bool b);
    static void WriteInt(FILE *fm, int i);
    static void WriteDouble(FILE *fm, double d);
    static void WriteUnsignedLong(FILE *fm, unsigned long ul);
    static void
      WriteUnsignedCharArray(FILE *fm, unsigned char *buffer, size_t length);
    static int WriteCompactSegmentation(FILE *fm,
                                        int *segmentation,
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
                      double zfpBitBudget);
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
                    double zfpBitBudget);

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
    void *inputData_;
    void *outputData_;
    Triangulation *triangulation_;
    TopologicalSimplification topologicalSimplification;

    // Parameters
    int compressionType_;
    bool zfpOnly_;
    int sqMethodInt_;
    std::string sqMethod_;
    int dataScalarType_;
    int dataExtent_[6];
    double dataSpacing_[3];
    double dataOrigin_[3];
    double tolerance_;
    double maximumError_;
    double zfpBitBudget_;
    bool dontSubdivide_;
    bool useTopologicalSimplification_;

    // Persistence compression.
    std::vector<int> segmentation_;
    std::vector<std::tuple<double, int>> mapping_;
    std::vector<std::tuple<int, double, int>> criticalConstraints_;

    // IO.
    int nbVertices;
    int nbSegments;
    int rawFileLength;
    std::vector<double> decompressedData_;
    std::vector<int> decompressedOffsets_;
    std::vector<int> compressedOffsets_;
    int vertexNumberRead_;
    char *fileName;
  };

  // End namespace ttk.
} // namespace ttk

#include <OtherCompression.h>
#include <PersistenceDiagramCompression.h>

template <class dataType>
int ttk::TopologicalCompression::execute(const double &tol) {
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Starting compression... " << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

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
                                             double zfpBitBudget) {
  // [->fp] Write metadata.
  WriteMetaData<double>(fp, compressionType, zfpOnly, sqMethod, dataType,
                        dataExtent, dataSpacing, dataOrigin, tolerance,
                        zfpBitBudget);

#ifdef TTK_ENABLE_ZLIB
  WriteBool(fp, true);
#else
  WriteBool(fp, false);
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
  std::stringstream str;
  str << fileName << ".temp";
  const std::string s = str.str();
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

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Topology successfully written to buffer."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

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
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Geometry successfully written to buffer."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  } else {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Geometry was not successfully written "
             "to buffer."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
    if(fflush(fp))
      fclose(fp);
    else
      fclose(fp);
    return -1;
  }

  // Check computed size vs read size.
  if(totalSize < rawFileLength) {
    std::stringstream msg;
    msg << "[TopologicalCompression] Invalid total size (" << totalSize
        << " vs " << rawFileLength << ")." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

#ifdef TTK_ENABLE_ZLIB
  // [fm->ff] Compress fm.
  uLong sourceLen = (uLong)rawFileLength;
  Bytef *source = reinterpret_cast<unsigned char *>(buf);
  uLongf destLen = compressBound(sourceLen);
  std::vector<Bytef> ddest(destLen);
  Bytef *dest = ddest.data();
  CompressWithZlib(false, dest, &destLen, source, sourceLen);
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Data successfully compressed."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  // [fm->fp] Copy fm to fp.
  WriteUnsignedLong(fp, destLen); // Compressed size...
  WriteUnsignedLong(fp, sourceLen);
  WriteUnsignedCharArray(fp, dest, destLen);
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Data successfully written to filesystem."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }
#else
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] ZLIB not found, writing raw file."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }
  unsigned char *source = reinterpret_cast<unsigned char *>(buf);
  unsigned long sourceLen = (unsigned long)rawFileLength;
  unsigned long destLen = (unsigned long)rawFileLength;

  WriteUnsignedLong(fp, destLen); // Compressed size...
  WriteUnsignedLong(fp, sourceLen);
  WriteUnsignedCharArray(fp, source, destLen);
#endif

  if(fflush(fp))
    fclose(fp);
  else
    fclose(fp);

  return ret;
}

template <typename T>
int ttk::TopologicalCompression::WriteMetaData(FILE *fp,
                                               int compressionType,
                                               bool zfpOnly,
                                               const char *sqMethod,
                                               int dataType,
                                               int *dataExtent,
                                               double *dataSpacing,
                                               double *dataOrigin,
                                               double tolerance,
                                               double zfpBitBudget) {
  // -2. Persistence, or Other
  WriteInt(fp, compressionType);

  // -1. zfpOnly
  WriteBool(fp, zfpOnly);

  // 0. SQ type
  const char *sq = sqMethod;
  int sqType = (strcmp(sq, "") == 0)
                 ? 0
                 : (strcmp(sq, "r") == 0 || strcmp(sq, "R") == 0)
                     ? 1
                     : (strcmp(sq, "d") == 0 || strcmp(sq, "D") == 0) ? 2 : 3;

  WriteInt(fp, sqType);

  // 1. DataType
  WriteInt(fp, dataType);

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    WriteInt(fp, dataExtent[i]);

  for(int i = 0; i < 3; ++i)
    WriteDouble(fp, dataSpacing[i]);

  for(int i = 0; i < 3; ++i)
    WriteDouble(fp, dataOrigin[i]);

  // 4. Tolerance
  WriteDouble(fp, tolerance);

  // 5. ZFP ratio
  WriteDouble(fp, zfpBitBudget);

  {
    std::stringstream msg;
    msg << "[ttkCompressionWriter] Metadata successfully written." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  return 0;
}

template <typename T>
int ttk::TopologicalCompression::ReadFromFile(FILE *fp) {
  // [fp->] Read headers.
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Successfully read metadata." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  if(zfpOnly_ && (zfpBitBudget_ > 64 || zfpBitBudget_ < 1)) {
    std::stringstream msg;
    msg << "[TopologicalCompression] Wrong ZFP bit budget for ZFP-only use."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    return -4;
  }

  bool useZlib = ReadBool(fp);
  unsigned char *dest;
  std::vector<unsigned char> ddest;
  unsigned long sourceLen;
  unsigned long destLen;

#ifdef TTK_ENABLE_ZLIB
  if(useZlib) {
    // [fp->ff] Read compressed data.
    uLongf sl = ReadUnsignedLong(fp); // Compressed size...
    uLongf dl = ReadUnsignedLong(fp); // Uncompressed size...

    sourceLen = (uLongf)sl;
    destLen = dl;
    std::vector<Bytef> ssource(sl);
    Bytef *source = ssource.data();
    ReadUnsignedCharArray(fp, source, sl);
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully read compressed data."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }

    // [ff->fm] Decompress data.
    ddest.resize(destLen);
    dest = ddest.data();
    CompressWithZlib(true, dest, &destLen, source, sourceLen);
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully uncompressed data."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  } else {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] File was not compressed with ZLIB."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }

    ReadUnsignedLong(fp); // Compressed size...
    unsigned long dl = ReadUnsignedLong(fp); // Uncompressed size...

    destLen = dl;
    ddest.resize(destLen);
    dest = ddest.data();
    ReadUnsignedCharArray(fp, dest, destLen);
  }
#else
  if(useZlib) {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] File compressed but ZLIB not installed! "
             "Aborting."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
    return -4;
  } else {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] ZLIB not installed, but file was not "
             "compressed anyways."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }

    ReadUnsignedLong(fp); // Compressed size...
    unsigned long dl = ReadUnsignedLong(fp); // Uncompressed size...

    destLen = dl;
    ddest.resize(destLen);
    dest = ddest.data();
    ReadUnsignedCharArray(fp, dest, destLen);
  }
#endif

  // [fm->] Read data.
  char *buf = reinterpret_cast<char *>(dest);

  //#ifndef _MSC_VER
  // FILE *fm = fmemopen(buf, destLen, "r+");
  //#else
  std::stringstream str;
  str << fileName << ".temp";
  const std::string s = str.str();
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

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Successfully read topology." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

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
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully read geometry."
          << std::endl;
      msg << "[TopologicalCompression] Successfully read file." << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  } else {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Failed to write (possibly ZFP)!"
          << std::endl;
      msg << "[TopologicalCompression] File may be corrupted!" << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  }

  return 0;
}

template <typename T>
int ttk::TopologicalCompression::ReadMetaData(FILE *fm) {
  // -2. Compression type.
  compressionType_ = ReadInt(fm);

  // -1. ZFP only type.
  zfpOnly_ = ReadBool(fm);

  // 0. SQ type
  sqMethodInt_ = ReadInt(fm);

  // 1. DataType
  dataScalarType_ = ReadInt(fm);
  // DataScalarType = VTK_DOUBLE;

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    dataExtent_[i] = ReadInt(fm);

  for(int i = 0; i < 3; ++i)
    dataSpacing_[i] = ReadDouble(fm);

  for(int i = 0; i < 3; ++i)
    dataOrigin_[i] = ReadDouble(fm);

  // 4. Error tolerance (relative percentage)
  tolerance_ = ReadDouble(fm);

  // 5. Lossy compressor ratio
  zfpBitBudget_ = ReadDouble(fm);

  return 0;
}

#endif // TOPOLOGICALCOMPRESSION_H
