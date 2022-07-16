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
#include <cstring>
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
    int computeOther() const;

    template <typename dataType>
    int compressForOther(int vertexNumber,
                         const dataType *const inputData,
                         const SimplexId *const inputOffsets,
                         dataType *outputData,
                         const double &tol) const;

    // Getters and setters.
    inline void setCompressionType(int compressionType) {
      compressionType_ = compressionType;
    }
    inline void setSQ(const std::string &sqMethod) {
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

    inline int getNbVertices() const {
      return NbVertices;
    }
    inline int getNbSegments() const {
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

    inline int getCompressionType() const {
      return compressionType_;
    }
    inline int getSQMethod() const {
      return SQMethodInt;
    }
    inline int getDataScalarType() const {
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
    inline double getTolerance() const {
      return Tolerance;
    }
    inline double getZFPTolerance() const {
      return ZFPTolerance;
    }
    /**
     * @brief Switch from a relative (% of scalar field range) to an
     * absolute ZFP tolerance
     */
    inline void relToAbsZFPTolerance(const double zfpRelTol,
                                     const std::array<double, 2> &sfRange) {
      this->ZFPTolerance = zfpRelTol * (sfRange[1] - sfRange[0]) / 100.0;
    }
    inline bool getZFPOnly() const {
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
    }
    inline static bool cmp2(const std::tuple<double, int> &a,
                            const std::tuple<double, int> &b) {
      return std::get<1>(a) < std::get<1>(b);
    }

    template <typename T>
    T Read(FILE *fm) const {
      static_assert(std::is_same<T, uint8_t>() || std::is_same<T, int32_t>()
                      || std::is_same<T, uint64_t>()
                      || std::is_same<T, double>(),
                    "Function should have only those types");
      T ret;
      const auto status = std::fread(&ret, sizeof(T), 1, fm);
      if(status == 0) {
        this->printErr("Error reading " + std::string(typeid(T).name()) + "!");
      }
      return ret;
    }
    template <typename T>
    void ReadByteArray(FILE *fm, T *buffer, size_t length) const {
      const auto status = std::fread(buffer, sizeof(T), length, fm);
      if(status == 0) {
        this->printErr("Error reading " + std::string(typeid(T).name())
                       + "array!");
      }
    }
    template <typename T>
    void Write(FILE *fm, T data) const {
      static_assert(std::is_same<T, uint8_t>() || std::is_same<T, int32_t>()
                      || std::is_same<T, uint64_t>()
                      || std::is_same<T, double>(),
                    "Function should have only those types");
      auto status = std::fwrite(&data, sizeof(T), 1, fm);
      if(status == 0) {
        this->printErr("Error writing " + std::string(typeid(T).name()) + "!");
      }
    }
    template <typename T>
    void WriteByteArray(FILE *fm, const T *buffer, size_t length) const {
      const auto status = std::fwrite(buffer, sizeof(T), length, fm);
      if(status == 0) {
        this->printErr("Error writing " + std::string(typeid(T).name())
                       + "array!");
      }
    }

    int ReadCompactSegmentation(FILE *fm,
                                std::vector<int> &segmentation,
                                int &numberOfVertices,
                                int &numberOfSegments) const;
    int ReadPersistenceIndex(
      FILE *fm,
      std::vector<std::tuple<double, int>> &mappings,
      std::vector<std::tuple<double, int>> &mappingsSortedPerValue,
      std::vector<std::tuple<int, double, int>> &constraints,
      double &min,
      double &max,
      int &nbConstraints) const;

    int ReadMetaData(FILE *fm);
    template <typename triangulationType>
    int ReadFromFile(FILE *fm, const triangulationType &triangulation);

    int WriteCompactSegmentation(FILE *fm,
                                 const std::vector<int> &segmentation,
                                 int numberOfVertices,
                                 int numberOfSegments) const;
    int WritePersistenceIndex(
      FILE *fm,
      std::vector<std::tuple<double, int>> &mapping,
      std::vector<std::tuple<int, double, int>> &constraints) const;

    int WriteMetaData(FILE *fp,
                      int compressionType,
                      bool zfpOnly,
                      const char *sqMethod,
                      int dataType,
                      int *dataExtent,
                      double *dataSpacing,
                      double *dataOrigin,
                      double tolerance,
                      double zfpTolerance,
                      const std::string &dataArrayName);

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
                    double zfpTolerance,
                    const std::string &dataArrayName);

    template <typename dataType>
    void CropIntervals(
      std::vector<std::tuple<dataType, int>> &mappings,
      std::vector<std::tuple<dataType, int>> &mappingsSortedPerValue,
      double min,
      double max,
      int vertexNumber,
      double *array,
      std::vector<int> &Seg) const;

    // API management.

#ifdef TTK_ENABLE_ZFP
    int CompressWithZFP(FILE *file,
                        const bool decompress,
                        std::vector<double> &array,
                        const int nx,
                        const int ny,
                        const int nz,
                        const double zfpTolerance) const;
#endif

#ifdef TTK_ENABLE_ZLIB
    unsigned long GetZlibDestLen(const unsigned long sourceLen) const;
    void CompressWithZlib(bool decompress,
                          unsigned char *dest,
                          unsigned long &destLen,
                          const unsigned char *const source,
                          const unsigned long sourceLen) const;
#endif

  private:
    // Internal read/write.

    int ComputeTotalSizeForOther() const;

    int ComputeTotalSizeForPersistenceDiagram(
      std::vector<std::tuple<double, int>> &mapping,
      std::vector<std::tuple<int, double, int>> &criticalConstraints,
      bool zfpOnly,
      int nbSegments,
      int nbVertices,
      double zfpTolerance) const;

    int ReadPersistenceTopology(FILE *fm);
    int ReadOtherTopology(FILE *fm) const;
    template <typename triangulationType>
    int ReadPersistenceGeometry(FILE *fm,
                                const triangulationType &triangulation);
    int ReadOtherGeometry(FILE *fm) const;

    int WritePersistenceTopology(FILE *fm);
    int WriteOtherTopology(FILE *fm) const;
    int WritePersistenceGeometry(FILE *fm,
                                 int *dataExtent,
                                 bool zfpOnly,
                                 double zfpTolerance,
                                 double *toCompress);
    int WriteOtherGeometry(FILE *fm) const;

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
    double ZFPTolerance{50};
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
    const unsigned long formatVersion_{2};
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

template <typename triangulationType>
int ttk::TopologicalCompression::ReadFromFile(
  FILE *fp, const triangulationType &triangulation) {
  // [fp->] Read headers.

  this->printMsg("Successfully read metadata.");

  if(ZFPOnly && ZFPTolerance < 0.0) {
    this->printMsg("Wrong ZFP absolute error tolerance for ZFP-only use.");
    return -4;
  }

  bool useZlib = Read<uint8_t>(fp);
  unsigned char *dest;
  std::vector<unsigned char> ddest;
  unsigned long destLen;

#ifdef TTK_ENABLE_ZLIB
  if(useZlib) {
    // [fp->ff] Read compressed data.
    auto sl = Read<uint64_t>(fp); // Compressed size...
    auto dl = Read<uint64_t>(fp); // Uncompressed size...

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

    Read<uint64_t>(fp); // Compressed size...
    const auto dl = Read<uint64_t>(fp); // Uncompressed size...

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

    Read<uint64_t>(fp); // Compressed size...
    const auto dl = Read<uint64_t>(fp); // Uncompressed size...

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
      ReadPersistenceTopology(fm);
    else if(compressionType_ == (int)ttk::CompressionType::Other)
      ReadOtherTopology(fm);
  }

  this->printMsg("Successfully read topology.");

  // Get altered geometry.
  // Rebuild topologically consistent geometry.
  int status = 0;
  if(compressionType_ == (int)ttk::CompressionType::PersistenceDiagram)
    status = ReadPersistenceGeometry(fm, triangulation);
  else if(compressionType_ == (int)ttk::CompressionType::Other)
    status = ReadOtherGeometry(fm);

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
