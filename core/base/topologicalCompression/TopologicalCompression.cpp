#include "TopologicalCompression.h"

// General.
ttk::TopologicalCompression::TopologicalCompression() {
  this->setDebugMsgPrefix("TopologicalCompression");
}

// Dependencies.

#ifdef TTK_ENABLE_ZFP

#include <zfp.h>

int ttk::TopologicalCompression::CompressWithZFP(
  FILE *file,
  const bool decompress,
  std::vector<double> &array,
  const int nx,
  const int ny,
  const int nz,
  const double zfpTolerance) const {

  int n1 = 0, n2 = 0;
  bool is2D = nx == 1 || ny == 1 || nz == 1;
  if(is2D) {
    if(nx + ny == 2 || ny + nz == 2 || nx + nz == 2) {
      this->printErr("One-dimensional arrays not supported.");
      return 0;
    }

    n1 = nx != 1 ? nx : ny;
    n2 = nx != 1 && ny != 1 ? ny : nz;
  }

  // array scalar type
  const auto type = zfp_type_double;

  // array meta data
  zfp_field *field;
  if(is2D) {
    field = zfp_field_2d(array.data(), type, n1, n2);
  } else {
    field = zfp_field_3d(array.data(), type, nx, ny, nz);
  }

  // compressed stream
  auto *zfp = zfp_stream_open(nullptr);

  // set compression mode and parameters via one of three functions (c.f.
  // https://zfp.readthedocs.io/en/release0.5.5/modes.html)
  // use fixed-accuracy mode
  zfp_stream_set_accuracy(zfp, zfpTolerance);

  // allocate buffer for compressed data
  const auto bufsize = zfp_stream_maximum_size(zfp, field);
  std::vector<unsigned char> buffer(bufsize);

  // associate bit stream with allocated buffer
  auto *stream = stream_open(buffer.data(), bufsize);
  zfp_stream_set_bit_stream(zfp, stream);
  zfp_stream_rewind(zfp);

  // return value: 0 = success
  int status = 0;

  // byte size of compressed stream
  size_t zfpsize{};

  // compress or decompress entire array
  if(decompress) {
    // read compressed stream and decompress array
    zfpsize += fread(buffer.data(), 1, bufsize, file);

    // read the ZFP header (from v2)
    const auto res = zfp_read_header(zfp, field, ZFP_HEADER_FULL);
    if(res == 0) {
      this->printErr("Could not read ZFP header");
      status = 1;
    }

    if(!zfp_decompress(zfp, field)) {
      this->printErr("Decompression failed");
      status = 1;
    }
  } else {
    // write the ZFP header
    const auto res = zfp_write_header(zfp, field, ZFP_HEADER_FULL);
    if(res == 0) {
      this->printErr("Could not write ZFP header");
      status = 1;
    }

    // compress array and output compressed stream
    zfpsize += zfp_compress(zfp, field);
    if(!zfpsize) {
      this->printErr("Compression failed");
      status = 1;
    } else
      fwrite(buffer.data(), 1, zfpsize, file);
  }

  // clean up
  zfp_field_free(field);
  zfp_stream_close(zfp);
  stream_close(stream);

  if(status != 0) {
    this->printErr("Encountered a problem with ZFP.");
  }

  return (int)zfpsize;
}

#endif // TTK_ENABLE_ZFP

#ifdef TTK_ENABLE_ZLIB

#include <zlib.h>

unsigned long ttk::TopologicalCompression::GetZlibDestLen(
  const unsigned long sourceLen) const {
  return compressBound(sourceLen);
}

void ttk::TopologicalCompression::CompressWithZlib(
  bool decompress,
  unsigned char *dest,
  unsigned long &destLen,
  const unsigned char *const source,
  const unsigned long sourceLen) const {

  if(decompress)
    uncompress(dest, &destLen, source, sourceLen);
  else
    compress(dest, &destLen, source, sourceLen);
}

#endif // TTK_ENABLE_ZLIB

unsigned int ttk::TopologicalCompression::log2(int val) {
  if(val == 0)
    return UINT_MAX;
  if(val == 1)
    return 0;
  unsigned int ret = 0;
  while(val > 1) {
    val >>= 1;
    ret++;
  }
  return ret;
}

// IO.

int ttk::TopologicalCompression::ReadCompactSegmentation(
  FILE *fm,
  std::vector<int> &segmentation,
  int &numberOfVertices,
  int &numberOfSegments) const {

  int numberOfBytesRead = 0;

  numberOfBytesRead += sizeof(int);
  numberOfVertices = Read<int32_t>(fm);

  numberOfBytesRead += sizeof(int);
  numberOfSegments = Read<int32_t>(fm);

  unsigned int numberOfBitsPerSegment = log2(numberOfSegments) + 1;

#ifndef TTK_ENABLE_KAMIKAZE
  // avoid left shift with negative operand
  if(numberOfBitsPerSegment == 0) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // [MEDIUM] TODO: support long int
  if(numberOfBitsPerSegment > 32)
    return -3;

  // Decode
  int currentCell = 0;
  int offset = 0;
  int maskerRank = 0;
  int oldCompressedInt = 0;

  while(currentCell < numberOfVertices) {

    int compressedInt;
    numberOfBytesRead += sizeof(int);
    compressedInt = Read<int32_t>(fm);

    while(offset + numberOfBitsPerSegment <= 32) {

      // Size of segment < size of compressed int
      int currentSegment = compressedInt;

      if(maskerRank == 0) {
        currentSegment <<= (32 - offset - numberOfBitsPerSegment);

        if(currentSegment < 0) {
          currentSegment &= 2147483647;
          currentSegment >>= (32 - numberOfBitsPerSegment);
          currentSegment |= (1 << (numberOfBitsPerSegment - 1));
        } else {
          currentSegment >>= (32 - numberOfBitsPerSegment);
        }

        offset += numberOfBitsPerSegment;
      }

      // Got overlapping mask.
      else {
        currentSegment <<= (32 - offset);
        if(currentSegment < 0) {
          currentSegment &= 2147483647;
          currentSegment >>= (32 - numberOfBitsPerSegment);
          currentSegment |= (1 << (numberOfBitsPerSegment - 1));
        } else {
          currentSegment >>= (32 - numberOfBitsPerSegment);
        }

        int nextSegment = currentSegment;

        if(oldCompressedInt < 0) {
          oldCompressedInt &= 2147483647;
          oldCompressedInt >>= maskerRank;
          oldCompressedInt |= (1 << (32 - maskerRank - 1));
        } else {
          oldCompressedInt >>= maskerRank;
        }

        currentSegment = nextSegment | oldCompressedInt;
        maskerRank = 0;
      }

      segmentation.push_back(currentSegment);
      currentCell++; // Go to next cell.
    }

    // Overlapping mask into next container.
    {
      oldCompressedInt = compressedInt;
      if(offset == 32) {
        maskerRank = 0;
      } else {
        maskerRank = offset;
        offset += numberOfBitsPerSegment;
      }
      offset %= 32;
    }
  }

  return numberOfBytesRead;
}

// Returns number of bytes written.
int ttk::TopologicalCompression::WriteCompactSegmentation(
  FILE *fm,
  const std::vector<int> &segmentation,
  int numberOfVertices,
  int numberOfSegments) const {

  int numberOfBytesWritten = 0;

  // Compute number of bits per segment
  // (can be deduced at read-time from numberOfSegments)
  unsigned int numberOfBitsPerSegment = log2(numberOfSegments) + 1;

  // [MEDIUM] TODO: support long int
  if(numberOfBitsPerSegment > 32)
    return -3;

  // Encode
  int currentCell = 0;
  int offset = 0;
  int maskerRank = 0;
  while(currentCell < numberOfVertices) {

    // Create next container.
    int compressedInt = 0;

    // Write regularly all segments until one overlaps
    // two containers.
    while(offset + numberOfBitsPerSegment <= 32) {

      // out-of-bounds here if segmentation.size() == numberOfVertices
      // (segmentation allocated in compressForPersistenceDiagram)
      int currentSegment = segmentation[currentCell];

      // If applicable, fill last part of current segment.
      if(maskerRank != 0) {
        // Always a positive int.
        currentSegment = currentSegment >> (numberOfBitsPerSegment - offset);
        // (numberOfBitsPerSegment - maskerRank - 1);
        maskerRank = 0;
        compressedInt |= currentSegment;

      } else {
        int cursor = (currentSegment << offset); // 0es after <<
        compressedInt |= cursor;
        offset += numberOfBitsPerSegment;
      }

      currentCell++; // Next cell.
    }

    // Test for overflow filling last part of current container.
    if(currentCell >= numberOfVertices) {
      numberOfBytesWritten += sizeof(int);
      Write<int32_t>(fm, compressedInt);
      break;
    }

    // Write current segment into last part of current container,
    // to be continued into next container.
    {
      int currentSegment = segmentation[currentCell];

      if(offset == 32) {
        // currentCell++;
      } else {
        int cursor = (currentSegment << offset);
        compressedInt = compressedInt | cursor;

        maskerRank = 32 - offset;
        offset += numberOfBitsPerSegment;
      }

      offset %= 32;
    }

    // Dump current container.
    numberOfBytesWritten += sizeof(int);
    Write<int32_t>(fm, compressedInt);
  }

  return numberOfBytesWritten;
}

int ttk::TopologicalCompression::ReadPersistenceIndex(
  FILE *fm,
  std::vector<std::tuple<double, int>> &mappings,
  std::vector<std::tuple<double, int>> &mappingsSortedPerValue,
  std::vector<std::tuple<int, double, int>> &constraints,
  double &min,
  double &max,
  int &nbConstraints) const {

  int numberOfBytesRead = 0;

  // 1.a. Read mapping.
  int mappingSize;
  numberOfBytesRead += sizeof(int);
  mappingSize = Read<int32_t>(fm);

  for(int i = 0; i < mappingSize; ++i) {
    int idv;
    numberOfBytesRead += sizeof(int);
    idv = Read<int32_t>(fm);

    double value;
    numberOfBytesRead += sizeof(double);
    value = Read<double>(fm);

    mappings.push_back(std::make_tuple(value, idv));
    mappingsSortedPerValue.push_back(std::make_tuple(value, idv));
  }

  // Sort mapping.
  std::sort(mappings.begin(), mappings.end(), cmp);
  std::sort(mappingsSortedPerValue.begin(), mappingsSortedPerValue.end(), cmp2);

  // 1.b. Read constraints.
  numberOfBytesRead += sizeof(int);
  nbConstraints = Read<int32_t>(fm);

  for(int i = 0; i < nbConstraints; ++i) {
    int idVertex;
    double value;
    int vertexType;

    numberOfBytesRead += sizeof(int);
    idVertex = Read<int32_t>(fm);

    numberOfBytesRead += sizeof(double);
    value = Read<double>(fm);

    numberOfBytesRead += sizeof(int);
    vertexType = Read<int32_t>(fm);

    if(i == 0) {
      min = value;
      max = value;
    }

    if(value < min)
      min = value;
    if(value > max)
      max = value;

    constraints.push_back(std::make_tuple(idVertex, value, vertexType));
  }

  return numberOfBytesRead;
}

int ttk::TopologicalCompression::WritePersistenceIndex(
  FILE *fm,
  std::vector<std::tuple<double, int>> &mapping,
  std::vector<std::tuple<int, double, int>> &constraints) const {

  int numberOfBytesWritten = 0;

  // Size.
  auto mappingSize = (int)mapping.size();
  numberOfBytesWritten += sizeof(int);
  Write<int32_t>(fm, mappingSize);

  // Segmentation values for each particular index.
  for(int i = 0; i < mappingSize; ++i) {
    std::tuple<double, int> t = mapping[i];
    int idv = std::get<1>(t);
    numberOfBytesWritten += sizeof(int);
    Write<int32_t>(fm, idv);

    auto value = std::get<0>(t);
    numberOfBytesWritten += sizeof(double);
    Write<double>(fm, value);
  }

  auto nbConstraints = (int)constraints.size();
  numberOfBytesWritten += sizeof(int);
  Write<int32_t>(fm, nbConstraints);

  for(int i = 0; i < nbConstraints; ++i) {
    std::tuple<int, double, int> t = constraints[i];
    int idVertex = std::get<0>(t);
    auto value = std::get<1>(t);
    int vertexType = std::get<2>(t);

    numberOfBytesWritten += sizeof(int);
    Write<int32_t>(fm, idVertex);

    numberOfBytesWritten += sizeof(double);
    Write<double>(fm, value);

    numberOfBytesWritten += sizeof(int);
    Write<int32_t>(fm, vertexType);
  }

  return numberOfBytesWritten;
}

int ttk::TopologicalCompression::ComputeTotalSizeForPersistenceDiagram(
  std::vector<std::tuple<double, int>> &mapping,
  std::vector<std::tuple<int, double, int>> &criticalConstraints,
  bool zfpOnly,
  int nSegments,
  int nVertices,
  double zfpTolerance) const {

  int totalSize = 0;

  if(!zfpOnly) {
    // Topological segments.
    int numberOfBitsPerSegment = log2(nSegments) + 1;
    double nbCharPerSegment = (double)numberOfBitsPerSegment / 8.0;
    totalSize += (sizeof(int) * 2 + std::ceil(nbCharPerSegment * nVertices));

    // Geometrical mapping.
    auto mappingSize = (int)mapping.size();
    auto constraintsSize = (int)criticalConstraints.size();
    totalSize += (mappingSize) * (sizeof(int) + sizeof(double)) + sizeof(int);
    totalSize
      += (constraintsSize) * (2 * sizeof(int) + sizeof(double)) + sizeof(int);
  }

  // conservative estimate of the ZFP buffer (no compression at all...)
  totalSize += zfpTolerance > 0.0 ? nVertices * sizeof(double) : 0 + 2;

  return totalSize;
}

int ttk::TopologicalCompression::WritePersistenceTopology(FILE *fm) {
  int numberOfBytesWritten = 0;

  int numberOfVertices = getNbVertices();
  int numberOfSegments = getNbSegments();

  // Test arguments.
  if(numberOfSegments < 1)
    return -1;

  numberOfBytesWritten += sizeof(int);
  Write<int32_t>(fm, numberOfVertices);

  numberOfBytesWritten += sizeof(int);
  Write<int32_t>(fm, numberOfSegments);

  numberOfBytesWritten += WriteCompactSegmentation(
    fm, getSegmentation(), numberOfVertices, numberOfSegments);

  rawFileLength += numberOfBytesWritten;

  return 0;
}

int ttk::TopologicalCompression::WritePersistenceGeometry(FILE *fm,
                                                          int *dataExtent,
                                                          bool zfpOnly,
                                                          double zfpTolerance,
                                                          double *toCompress) {
  int numberOfBytesWritten = 0;

  if(!zfpOnly) {
    // 1. Write segmentation map.
    // 2. Write critical constraints.
    numberOfBytesWritten
      += WritePersistenceIndex(fm, mapping_, criticalConstraints_);
  }

  this->printMsg("Wrote raw geometry.");

  if(zfpTolerance >= 0.0) {
#ifdef TTK_ENABLE_ZFP
    // (1. or 3.) Write zfp-compressed array.
    int nx = 1 + dataExtent[1] - dataExtent[0];
    int ny = 1 + dataExtent[3] - dataExtent[2];
    int nz = 1 + dataExtent[5] - dataExtent[4];

    std::vector<double> dataVector(toCompress, toCompress + (nx * ny * nz));
    numberOfBytesWritten
      += CompressWithZFP(fm, false, dataVector, nx, ny, nz, zfpTolerance);

#else
    TTK_FORCE_USE(dataExtent);
    TTK_FORCE_USE(toCompress);

    this->printErr("Attempted to write with ZFP but ZFP is not installed.");
    return -5;
#endif
  }

  rawFileLength += numberOfBytesWritten;

  return 0;
}

int ttk::TopologicalCompression::ReadPersistenceTopology(FILE *fm) {
  int numberOfSegments;
  int numberOfVertices;

  int numberOfBytesRead = ReadCompactSegmentation(
    fm, segmentation_, numberOfVertices, numberOfSegments);

  this->rawFileLength += numberOfBytesRead;

  return 0;
}

///////////////////////////////
// Other compression methods //
///////////////////////////////

int ttk::TopologicalCompression::ComputeTotalSizeForOther() const {
  // Should return the number of bytes to be written on the output file
  // sizeof(char) = 1 (byte)
  // use sizeof(int), sizeof(double) to get the number of bytes of
  // the matching structures.
  return 0;
}

int ttk::TopologicalCompression::computeOther() const {
  // Code me
  return 0;
}

int ttk::TopologicalCompression::WriteOtherTopology(
  FILE *ttkNotUsed(fm)) const {
  this->printWrn("Writing Other index / topology.");
  // Code me
  return 0;
}

int ttk::TopologicalCompression::WriteOtherGeometry(
  FILE *ttkNotUsed(fm)) const {
  this->printWrn("Writing Other buffer / geometry.");
  // Code me
  return 0;
}

int ttk::TopologicalCompression::ReadOtherTopology(FILE *ttkNotUsed(fm)) const {
  this->printWrn("Reading Other index / topology.");
  // Code me
  return 0;
}

int ttk::TopologicalCompression::ReadOtherGeometry(FILE *ttkNotUsed(fm)) const {
  this->printWrn("Reading Other buffer / geometry.");
  // Code me
  return 0;
}

/////////////////////////////
// Read/Write File methods //
/////////////////////////////

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
                                             double zfpTolerance,
                                             const std::string &dataArrayName) {
  // [->fp] Write metadata.
  WriteMetaData(fp, compressionType, zfpOnly, sqMethod, dataType, dataExtent,
                dataSpacing, dataOrigin, tolerance, zfpTolerance,
                dataArrayName);

#ifdef TTK_ENABLE_ZLIB
  Write<uint8_t>(fp, true);
#else
  Write<uint8_t>(fp, false);
#endif

  bool usePersistence
    = compressionType == (int)ttk::CompressionType::PersistenceDiagram;
  bool useOther = compressionType == (int)ttk::CompressionType::Other;

  int numberOfVertices = 1;
  for(int i = 0; i < 3; ++i)
    numberOfVertices *= (1 + dataExtent[2 * i + 1] - dataExtent[2 * i]);
  NbVertices = numberOfVertices;

  int totalSize = usePersistence ? ComputeTotalSizeForPersistenceDiagram(
                    getMapping(), getCriticalConstraints(), zfpOnly,
                    getNbSegments(), getNbVertices(), zfpTolerance)
                  : useOther ? ComputeTotalSizeForOther()
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
      WritePersistenceTopology(fm);
    else if(useOther)
      WriteOtherTopology(fm);
  }

  this->printMsg("Topology successfully written to buffer.");

  int status = 0;
  // [->fm] Write altered geometry.
  if(usePersistence)
    status
      = WritePersistenceGeometry(fm, dataExtent, zfpOnly, zfpTolerance, data);
  else if(useOther)
    status = WriteOtherGeometry(fm);

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
  Write<uint64_t>(fp, destLen); // Compressed size...
  Write<uint64_t>(fp, sourceLen);
  WriteByteArray(fp, ddest.data(), destLen);
  this->printMsg("Data successfully written to filesystem.");

#else
  this->printMsg("ZLIB not found, writing raw file.");
  unsigned char *source = reinterpret_cast<unsigned char *>(buf);
  unsigned long sourceLen = (unsigned long)rawFileLength;
  unsigned long destLen = (unsigned long)rawFileLength;

  Write<uint64_t>(fp, destLen); // Compressed size...
  Write<uint64_t>(fp, sourceLen);
  WriteByteArray(fp, source, destLen);
#endif

  fflush(fp);
  fclose(fp);

  return ret;
}

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
  double zfpTolerance,
  const std::string &dataArrayName) {

  // -4. Magic bytes
  WriteByteArray(fp, magicBytes_, std::strlen(magicBytes_));

  // -3. File format version
  Write<uint64_t>(fp, formatVersion_);

  // -2. Persistence, or Other
  Write<int32_t>(fp, compressionType);

  // -1. zfpOnly
  Write<uint8_t>(fp, zfpOnly);

  // 0. SQ type
  const char *sq = sqMethod;
  int sqType = (strcmp(sq, "") == 0)                            ? 0
               : (strcmp(sq, "r") == 0 || strcmp(sq, "R") == 0) ? 1
               : (strcmp(sq, "d") == 0 || strcmp(sq, "D") == 0) ? 2
                                                                : 3;

  Write<int32_t>(fp, sqType);

  // 1. DataType
  Write<int32_t>(fp, dataType);

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    Write<int32_t>(fp, dataExtent[i]);

  for(int i = 0; i < 3; ++i)
    Write<double>(fp, dataSpacing[i]);

  for(int i = 0; i < 3; ++i)
    Write<double>(fp, dataOrigin[i]);

  // 4. Tolerance
  Write<double>(fp, tolerance);

  // 5. ZFP ratio
  Write<double>(fp, zfpTolerance);

  // 6. Length of array name
  // (explicit call to unsigned long variant for MSVC compatibility)
  Write<uint64_t>(fp, dataArrayName.size());

  // 7. Array name (as unsigned chars)
  WriteByteArray(fp, dataArrayName.c_str(), dataArrayName.size());

  this->printMsg("Metadata successfully written.");

  return 0;
}

int ttk::TopologicalCompression::ReadMetaData(FILE *fm) {

  // -4. Magic bytes
  const auto magicBytesLen{std::strlen(magicBytes_)};
  std::vector<char> mBytes(magicBytesLen + 1);
  mBytes[magicBytesLen] = '\0'; // NULL-termination
  ReadByteArray(fm, mBytes.data(), magicBytesLen);

  // To deal with pre-v1 file format (without scalar field array name)
  const bool hasMagicBytes = strcmp(mBytes.data(), magicBytes_) == 0;
  if(!hasMagicBytes) {
    this->printErr("Could not find magic bytes in input file!");
    return 1;
  }

  // -3. File format version
  const auto fileVersion = Read<uint64_t>(fm);
  if(fileVersion < this->formatVersion_) {
    this->printErr("Old format version detected (" + std::to_string(fileVersion)
                   + " vs. " + std::to_string(this->formatVersion_) + ").");
    this->printErr("Older formats are not supported!");
    return 1;
  } else if(fileVersion > this->formatVersion_) {
    this->printErr("Newer format version detected ("
                   + std::to_string(fileVersion) + " vs. "
                   + std::to_string(this->formatVersion_) + ").");
    this->printErr("Cannot read file with current TTK, try with to update.");
    return 1;
  }

  // -2. Compression type.
  compressionType_ = Read<int32_t>(fm);

  // -1. ZFP only type.
  ZFPOnly = Read<uint8_t>(fm);

  // 0. SQ type
  SQMethodInt = Read<int32_t>(fm);

  // 1. DataType
  dataScalarType_ = Read<int32_t>(fm);
  // DataScalarType = VTK_DOUBLE;

  // 2. Data extent, spacing, origin
  for(int i = 0; i < 6; ++i)
    dataExtent_[i] = Read<int32_t>(fm);

  for(int i = 0; i < 3; ++i)
    dataSpacing_[i] = Read<double>(fm);

  for(int i = 0; i < 3; ++i)
    dataOrigin_[i] = Read<double>(fm);

  // 4. Error tolerance (relative percentage)
  Tolerance = Read<double>(fm);

  // 5. Lossy compressor ratio
  ZFPTolerance = Read<double>(fm);

  // 6. Length of array name
  size_t dataArrayNameLength = Read<uint64_t>(fm);

  // 7. Array name (as unsigned chars)
  dataArrayName_.resize(dataArrayNameLength + 1);
  dataArrayName_[dataArrayNameLength] = '\0'; // NULL-termination
  ReadByteArray(fm, dataArrayName_.data(), dataArrayNameLength);

  return 0;
}
