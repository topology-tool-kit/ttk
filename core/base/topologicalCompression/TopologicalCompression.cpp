#include "TopologicalCompression.h"

// General.
ttk::TopologicalCompression::TopologicalCompression() {
  inputData_ = nullptr;
  outputData_ = nullptr;
  triangulation_ = nullptr;
  sqMethod_ = "";
  nbSegments = 0;
  nbVertices = 0;
  rawFileLength = 0;
}

ttk::TopologicalCompression::~TopologicalCompression() {
}

// Dependencies.

#ifdef TTK_ENABLE_ZFP

int ttk::TopologicalCompression::CompressWithZFP(FILE *file,
                                                 bool decompress,
                                                 std::vector<double> &array,
                                                 int nx,
                                                 int ny,
                                                 int nz,
                                                 double rate) {
  return ttk::TopologicalCompression::compressZFPInternal(
    array.data(), nx, ny, nz, rate, decompress, file);
}

int ttk::TopologicalCompression::compressZFPInternal(double *array,
                                                     int nx,
                                                     int ny,
                                                     int nz,
                                                     double rate,
                                                     bool decompress,
                                                     FILE *file) {
  int status = 0; // return value: 0 = success
  zfp_type type; // array scalar type
  zfp_field *field; // array meta data
  zfp_stream *zfp; // compressed stream
  size_t bufsize; // byte size of compressed buffer
  bitstream *stream; // bit stream to write to or read from
  size_t zfpsize; // byte size of compressed stream

  int n1 = 0;
  int n2 = 0;
  bool is2D = nx == 1 || ny == 1 || nz == 1;
  if(is2D) {
    if(nx + ny == 2 || ny + nz == 2 || nx + nz == 2) {
      fprintf(stderr, "One-dimensional arrays not supported.\n");
      return 0;
    }

    n1 = nx != 1 ? nx : ny;
    n2 = nx != 1 && ny != 1 ? ny : nz;
  }

  // allocate meta data for the 3D array a[nz][ny][nx]
  type = zfp_type_double;

  if(is2D) {
    field = zfp_field_2d(array, type, (unsigned int)n1, (unsigned int)n2);
  } else {
    field = zfp_field_3d(
      array, type, (unsigned int)nx, (unsigned int)ny, (unsigned int)nz);
  }

  // allocate meta data for a compressed stream
  zfp = zfp_stream_open(NULL);

  // set compression mode and parameters via one of three functions
  zfp_stream_set_rate(zfp, rate, type, 3, 0);
  //  zfp_stream_set_precision(zfp, precision);
  // zfp_stream_set_accuracy(zfp, tolerance);

  // allocate buffer for compressed data
  bufsize = zfp_stream_maximum_size(zfp, field);
  std::vector<unsigned char> buffer(bufsize);

  // associate bit stream with allocated buffer
  stream = stream_open(buffer.data(), bufsize);
  zfp_stream_set_bit_stream(zfp, stream);
  zfp_stream_rewind(zfp);

  // compress or decompress entire array
  if(decompress) {
    // read compressed stream and decompress array
    // zfpsize = fread(buffer, 1, bufsize, stdin);
    zfpsize = fread(buffer.data(), 1, bufsize, file);
    if(!zfp_decompress(zfp, field)) {
      fprintf(stderr, "decompression failed\n");
      status = 1;
    }
  } else {
    // compress array and output compressed stream
    zfpsize = zfp_compress(zfp, field);
    if(!zfpsize) {
      fprintf(stderr, "compression failed\n");
      status = 1;
    } else
      fwrite(buffer.data(), 1, zfpsize, file);
    // fwrite(buffer, 1, zfpsize, stdout);
  }

  // clean up
  zfp_field_free(field);
  zfp_stream_close(zfp);
  stream_close(stream);
  // free(array);

  if(status != 0) {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[TopologicalCompression] Encountered a problem with ZFP."
        << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::fatalMsg);
  }

  return (int)zfpsize;
}

#endif

#ifdef TTK_ENABLE_ZLIB

void ttk::TopologicalCompression::CompressWithZlib(bool decompress,
                                                   Bytef *dest,
                                                   uLongf *destLen,
                                                   const Bytef *source,
                                                   uLong sourceLen) {
  if(decompress)
    uncompress(dest, destLen, source, sourceLen);
  else
    compress(dest, destLen, source, sourceLen);
}

#endif

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

bool ttk::TopologicalCompression::ReadBool(FILE *fm) {
  bool b;
  int ret = (int)std::fread(&b, sizeof(bool), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error reading boolean!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
  return b;
}

void ttk::TopologicalCompression::WriteBool(FILE *fm, bool b) {
  int ret = (int)std::fwrite(&b, sizeof(bool), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error writing boolean!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

int ttk::TopologicalCompression::ReadInt(FILE *fm) {
  int i;
  int ret = (int)std::fread(&i, sizeof(int), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error reading integer!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
  return i;
}

void ttk::TopologicalCompression::WriteInt(FILE *fm, int i) {
  int ret = (int)std::fwrite(&i, sizeof(int), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error writing integer!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

double ttk::TopologicalCompression::ReadDouble(FILE *fm) {
  double d;
  int ret = (int)std::fread(&d, sizeof(double), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug dbg;
    msg << "[TopologicalCompression] Error reading double!" << std::endl;
    dbg.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
  return d;
}

void ttk::TopologicalCompression::WriteDouble(FILE *fm, double d) {
  int ret = (int)std::fwrite(&d, sizeof(double), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug dbg;
    msg << "[TopologicalCompression] Error writing double!" << std::endl;
    dbg.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

unsigned long ttk::TopologicalCompression::ReadUnsignedLong(FILE *fm) {
  unsigned long ul;
  int ret = (int)std::fread(&ul, sizeof(unsigned long), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error reading long!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
  return ul;
}

void ttk::TopologicalCompression::WriteUnsignedLong(FILE *fm,
                                                    unsigned long ul) {
  int ret = (int)std::fwrite(&ul, sizeof(unsigned long), 1, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error writing long!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

void ttk::TopologicalCompression::ReadUnsignedCharArray(FILE *fm,
                                                        unsigned char *buffer,
                                                        size_t length) {
  int ret = (int)std::fread(buffer, sizeof(unsigned char), length, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error reading char array!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

void ttk::TopologicalCompression::WriteUnsignedCharArray(FILE *fm,
                                                         unsigned char *buffer,
                                                         size_t length) {
  int ret = (int)std::fwrite(buffer, sizeof(unsigned char), length, fm);
  if(!ret) {
    std::stringstream msg;
    ttk::Debug d;
    msg << "[TopologicalCompression] Error writing char array!" << std::endl;
    d.dMsg(std::cerr, msg.str(), ttk::Debug::fatalMsg);
  }
}

int ttk::TopologicalCompression::ReadCompactSegmentation(
  FILE *fm,
  std::vector<int> &segmentation,
  int &numberOfVertices,
  int &numberOfSegments) {
  auto ReadInt = ttk::TopologicalCompression::ReadInt;

  int numberOfBytesRead = 0;

  numberOfBytesRead += sizeof(int);
  numberOfVertices = ReadInt(fm);

  numberOfBytesRead += sizeof(int);
  numberOfSegments = ReadInt(fm);

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
    compressedInt = ReadInt(fm);

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
  FILE *fm, int *segmentation, int numberOfVertices, int numberOfSegments) {
  auto WriteInt = ttk::TopologicalCompression::WriteInt;

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
      WriteInt(fm, compressedInt);
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
    WriteInt(fm, compressedInt);
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
  int &nbConstraints) {
  int numberOfBytesRead = 0;

  // 1.a. Read mapping.
  int mappingSize;
  numberOfBytesRead += sizeof(int);
  mappingSize = ReadInt(fm);

  for(int i = 0; i < mappingSize; ++i) {
    int idv;
    numberOfBytesRead += sizeof(int);
    idv = ReadInt(fm);

    double value;
    numberOfBytesRead += sizeof(double);
    value = ReadDouble(fm);

    mappings.push_back(std::make_tuple(value, idv));
    mappingsSortedPerValue.push_back(std::make_tuple(value, idv));
  }

  // Sort mapping.
  std::sort(mappings.begin(), mappings.end(), cmp);
  std::sort(mappingsSortedPerValue.begin(), mappingsSortedPerValue.end(), cmp2);

  // 1.b. Read constraints.
  numberOfBytesRead += sizeof(int);
  nbConstraints = ReadInt(fm);

  for(int i = 0; i < nbConstraints; ++i) {
    int idVertex;
    double value;
    int vertexType;

    numberOfBytesRead += sizeof(int);
    idVertex = ReadInt(fm);

    numberOfBytesRead += sizeof(double);
    value = ReadDouble(fm);

    numberOfBytesRead += sizeof(int);
    vertexType = ReadInt(fm);

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
  std::vector<std::tuple<int, double, int>> &constraints) {
  int numberOfBytesWritten = 0;

  // Size.
  auto mappingSize = (int)mapping.size();
  numberOfBytesWritten += sizeof(int);
  WriteInt(fm, mappingSize);

  // Segmentation values for each particular index.
  for(int i = 0; i < mappingSize; ++i) {
    std::tuple<double, int> t = mapping[i];
    int idv = std::get<1>(t);
    numberOfBytesWritten += sizeof(int);
    WriteInt(fm, idv);

    auto value = (double)std::get<0>(t);
    numberOfBytesWritten += sizeof(double);
    WriteDouble(fm, value);
  }

  auto nbConstraints = (int)constraints.size();
  numberOfBytesWritten += sizeof(int);
  WriteInt(fm, nbConstraints);

  for(int i = 0; i < nbConstraints; ++i) {
    std::tuple<int, double, int> t = constraints[i];
    int idVertex = std::get<0>(t);
    auto value = (double)std::get<1>(t);
    int vertexType = std::get<2>(t);

    numberOfBytesWritten += sizeof(int);
    WriteInt(fm, idVertex);

    numberOfBytesWritten += sizeof(double);
    WriteDouble(fm, value);

    numberOfBytesWritten += sizeof(int);
    WriteInt(fm, vertexType);
  }

  return numberOfBytesWritten;
}
