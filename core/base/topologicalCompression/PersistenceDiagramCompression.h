//
// Created by max on 24/05/18.
//

#ifndef TTK_PERSISTENCEDIAGRAMCOMPRESSION_H
#define TTK_PERSISTENCEDIAGRAMCOMPRESSION_H

template <typename dataType>
int ttk::TopologicalCompression::ComputeTotalSizeForPersistenceDiagram(
  std::vector<std::tuple<double, int>> &mapping,
  std::vector<std::tuple<int, double, int>> &criticalConstraints,
  bool zfpOnly,
  int nSegments,
  int nVertices,
  double zfpBitBudget) {
  using ttk::TopologicalCompression;

  int totalSize = 0;

  if(!(zfpOnly)) {
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

  totalSize += (zfpBitBudget <= 64 && zfpBitBudget > 0)
                 ? (nVertices * std::ceil(zfpBitBudget / 2.0))
                 : 0;
  totalSize += 2;
  // * 16 -> 32 bpd
  // -> * bitbudget / 2
  // -> * bitbudget / 2
  return totalSize;
}

template <typename dataType>
int ttk::TopologicalCompression::WritePersistenceTopology(FILE *fm) {
  int numberOfBytesWritten = 0;

  int numberOfVertices = getNbVertices();
  int numberOfSegments = getNbSegments();
  std::vector<int> segmentation = getSegmentation();

  // Test arguments.
  if(numberOfSegments < 1)
    return -1;

  numberOfBytesWritten += sizeof(int);
  WriteInt(fm, numberOfVertices);

  numberOfBytesWritten += sizeof(int);
  WriteInt(fm, numberOfSegments);

  numberOfBytesWritten += WriteCompactSegmentation(
    fm, segmentation.data(), numberOfVertices, numberOfSegments);

  rawFileLength += numberOfBytesWritten;

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::WritePersistenceGeometry(FILE *fm,
                                                          int *dataExtent,
                                                          bool zfpOnly,
                                                          double zfpBitBudget,
                                                          double *toCompress) {
  int numberOfBytesWritten = 0;

  if(!zfpOnly) {
    // 1. Write segmentation map.
    // 2. Write critical constraints.
    numberOfBytesWritten
      += WritePersistenceIndex(fm, mapping_, criticalConstraints_);
  }

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Wrote raw geometry." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  if(zfpBitBudget <= 64.0 && zfpBitBudget > 0) {
#ifdef TTK_ENABLE_ZFP
    // (1. or 3.) Write zfp-compressed array.
    int nx = 1 + dataExtent[1] - dataExtent[0];
    int ny = 1 + dataExtent[3] - dataExtent[2];
    int nz = 1 + dataExtent[5] - dataExtent[4];

    std::vector<double> dataVector(toCompress, toCompress + (nx * ny * nz));
    using ttk::TopologicalCompression;
    numberOfBytesWritten
      += CompressWithZFP(fm, false, dataVector, nx, ny, nz, zfpBitBudget);

#else
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Attempted to write with ZFP but ZFP is "
             "not installed."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
    return -5;
#endif
  }

  rawFileLength += numberOfBytesWritten;

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadPersistenceTopology(FILE *fm) {
  int numberOfSegments;
  int numberOfVertices;

  int numberOfBytesRead = ReadCompactSegmentation(
    fm, segmentation_, numberOfVertices, numberOfSegments);

  rawFileLength += numberOfBytesRead;

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadPersistenceGeometry(FILE *fm) {
  using ttk::TopologicalCompression;

  int sqMethod = sqMethodInt_;
  bool zfpOnly = zfpOnly_;
  double zfpBitBudget = zfpBitBudget_;
  int *dataExtent = dataExtent_;

  std::vector<std::tuple<double, int>> mappingsSortedPerValue;

  double min = 0;
  double max = 0;
  int nbConstraints = 0;

  int numberOfBytesRead = 0;
  if(!zfpOnly) {
    numberOfBytesRead
      += ReadPersistenceIndex(fm, mapping_, mappingsSortedPerValue,
                              criticalConstraints_, min, max, nbConstraints);
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully read geomap." << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  }

  // Prepare array reconstruction.
  int nx = 1 + dataExtent[1] - dataExtent[0];
  int ny = 1 + dataExtent[3] - dataExtent[2];
  int nz = 1 + dataExtent[5] - dataExtent[4];
  int vertexNumber = nx * ny * nz;

  decompressedData_.resize(vertexNumber);
  if(zfpBitBudget > 64.0 || zfpBitBudget < 1) {

    // 2.a. (2.) Affect values to points thanks to topology indices.
    for(int i = 0; i < vertexNumber; ++i) {
      int seg = segmentation_[i];
      auto end = mapping_.end();
      auto it = std::lower_bound(
        mapping_.begin(), mapping_.end(), std::make_tuple(0, seg), cmp);
      if(it != end) {
        std::tuple<double, int> tt = *it;
        double value = std::get<0>(tt);
        int sseg = std::get<1>(tt);
        if(seg != sseg) {
          std::stringstream msg;
          msg << "Decompression mismatch (" << seg << ", " << sseg << ")"
              << std::endl;
          dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
        }
        decompressedData_[i] = value;
      } else {
        {
          std::stringstream msg;
          msg << "Could not find " << seg << " index." << std::endl;
          dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
        }
        std::tuple<double, int> tt = *it;
        double value = std::get<0>(tt);
        decompressedData_[i] = value;
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully affected geomap."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }

  } else {
#ifdef TTK_ENABLE_ZFP
    // 2.b. (2.) Read with ZFP.
    using ttk::TopologicalCompression;
    numberOfBytesRead += zfpBitBudget * vertexNumber;
    CompressWithZFP(fm, true, decompressedData_, nx, ny, nz, zfpBitBudget);
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Successfully read with ZFP."
          << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
#else
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Attempted to read "
          << "a ZFP block but ZFP is not installed." << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
    return -5;
#endif
  }

  // No SQ.
  if(sqMethod == 0 || sqMethod == 3) {
    for(int i = 0; i < (int)criticalConstraints_.size(); ++i) {
      std::tuple<int, double, int> t = criticalConstraints_[i];
      int id = std::get<0>(t);
      double val = std::get<1>(t);
      decompressedData_[id] = val;
    }
  }

  // double tolerance = Tolerance;
  if(min == max) {
    std::stringstream msg;
    msg << "[TopologicalCompression] empty scalar field range." << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  } else {
    // tolerance *= (max - min);
  }

  if(sqMethod == 1 || sqMethod == 2)
    return 0;

  if(zfpOnly)
    return 0;

  // 2.b. (3.) Crop whatever doesn't fit in topological intervals.
  CropIntervals(mapping_, mappingsSortedPerValue, min, max, vertexNumber,
                decompressedData_.data(), segmentation_);
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Successfully cropped bad intervals."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  // 2.b. (4.) Apply topological simplification with min/max constraints
  PerformSimplification<double>(criticalConstraints_, nbConstraints,
                                vertexNumber, decompressedData_.data());
  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Successfully performed simplification."
        << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  rawFileLength += numberOfBytesRead;

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::PerformSimplification(
  const std::vector<std::tuple<int, double, int>> &constraints,
  int nbConstraints,
  int vertexNumber,
  double *array) {
  std::vector<int> inputOffsets(vertexNumber);
  std::vector<int> critConstraints(nbConstraints);
  std::vector<double> inArray(vertexNumber);
  // std::vector<int> *oo = new std::vector<int>(vertexNumber);
  decompressedOffsets_.resize(vertexNumber); // oo->data();
  int status = 0;

  // Offsets
  for(int i = 0; i < vertexNumber; ++i)
    inputOffsets[i] = i;

  // Triangulate.
  topologicalSimplification.setupTriangulation(triangulation_);

  // Preprocess simplification.
  // std::vector<int>* authorizedSaddles = new std::vector<int>();
  for(int i = 0; i < nbConstraints; ++i) {
    std::tuple<int, double, int> t = constraints[i];
    int id = std::get<0>(t);
    double val = std::get<1>(t);
    int type = std::get<2>(t);

    // if (type == 0)
    // authorizedSaddles->push_back(id);

    array[id] = val;

    // Smoothe neighborhood (along with offsets).
    SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(id);
    for(SimplexId j = 0; j < neighborNumber; ++j) {
      SimplexId neighbor;
      triangulation_->getVertexNeighbor(id, j, neighbor);

      if(type == 1) { // Local_maximum.
        if(array[neighbor] > val)
          array[neighbor] = val;
        if(array[neighbor] == val
           && inputOffsets[neighbor] > inputOffsets[id]) {
          int tmp = inputOffsets[id];
          inputOffsets[id] = inputOffsets[neighbor];
          inputOffsets[neighbor] = tmp;
        }
      } else if(type == -1) { // Local_minimum.
        if(array[neighbor] < val)
          array[neighbor] = val;
        if(array[neighbor] == val
           && inputOffsets[neighbor] < inputOffsets[id]) {
          int tmp = inputOffsets[id];
          inputOffsets[id] = inputOffsets[neighbor];
          inputOffsets[neighbor] = tmp;
        }
      } else if(type == 0) { // Saddle
      }
    }

    critConstraints[i] = id;
  }

  for(int i = 0; i < vertexNumber; ++i)
    inArray[i] = array[i];
  for(int i = 0; i < vertexNumber; ++i)
    decompressedOffsets_[i] = 0;

  topologicalSimplification.setInputScalarFieldPointer(inArray.data());
  topologicalSimplification.setOutputScalarFieldPointer(array);
  topologicalSimplification.setInputOffsetScalarFieldPointer(
    inputOffsets.data());
  topologicalSimplification.setOutputOffsetScalarFieldPointer(
    decompressedOffsets_.data());
  topologicalSimplification.setVertexIdentifierScalarFieldPointer(
    critConstraints.data());
  topologicalSimplification.setConstraintNumber(nbConstraints);
  status = topologicalSimplification.execute<double, int>();

  return status;
}

template <typename dataType>
void ttk::TopologicalCompression::CropIntervals(
  std::vector<std::tuple<dataType, int>> &mappings,
  std::vector<std::tuple<dataType, int>> &mappingsSortedPerValue,
  double min,
  double max,
  int vertexNumber,
  double *array,
  std::vector<int> &segmentation) {
  int numberOfMisses = 0;
  for(int i = 0; i < vertexNumber; ++i) {
    int seg = segmentation[i];
    auto end = mappings.end();
    auto it = lower_bound(
      mappings.begin(), mappings.end(), std::make_tuple(0, seg), cmp);
    if(it != end) {
      std::tuple<dataType, int> tt = *it;
      double value = std::get<0>(tt);
      int sseg = std::get<1>(tt);
      if(seg != sseg) {
        ttk::Debug d;
        std::stringstream msg;
        msg << "[TopologicalCompression] Decompression mismatch." << std::endl;
        d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
      }

      auto it2 = lower_bound(mappingsSortedPerValue.begin(),
                             mappingsSortedPerValue.end(),
                             std::make_tuple(value, 0), cmp2);

      if(it2 != mappingsSortedPerValue.end()
         && it2 != mappingsSortedPerValue.begin()
         && (it2 + 1) != mappingsSortedPerValue.end()) {
        dataType vv0 = std::get<0>(*(it2 - 1));
        // double vv1 = std::get<0>(*(it2));
        dataType vv2 = std::get<0>(*(++it2));
        // double m0 = (vv0 + vv1) / 2;
        // double m1 = (vv1 + vv2) / 2;

        if(array[i] < vv0) {
          numberOfMisses++;
          array[i] = vv0;
        } else if(array[i] > vv2) {
          numberOfMisses++;
          array[i] = vv2;
        }
      } else {
        if((it2 == mappingsSortedPerValue.end()
            || (it2 + 1) == mappingsSortedPerValue.end())
           && array[i] > max) {
          numberOfMisses++;
          array[i] = max;
        } else if((it2 == mappingsSortedPerValue.begin()
                   || (it2 - 1) == mappingsSortedPerValue.begin())
                  && array[i] < min) {
          numberOfMisses++;
          array[i] = min;
        }
      }
    } else {
      ttk::Debug d;
      std::stringstream msg;
      msg << "[TopologicalCompression] Error looking for topo index."
          << std::endl;
      d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  }

  if(numberOfMisses > 0) {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[TopologicalCompression] Missed " << numberOfMisses << " values."
        << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }
}

template <typename dataType>
int ttk::TopologicalCompression::computePersistencePairs(
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &JTPairs,
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &STPairs,
  dataType *inputScalars_,
  SimplexId *inputOffsets) {
  // Compute offsets
  const SimplexId numberOfVertices = triangulation_->getNumberOfVertices();
  std::vector<SimplexId> voffsets((unsigned long)numberOfVertices);
  std::copy(inputOffsets, inputOffsets + numberOfVertices, voffsets.begin());

  // Get contour tree
  ftm::FTMTreePP ftmTreePP;
  ftmTreePP.setupTriangulation(triangulation_, false);
  ftmTreePP.setVertexScalars(inputScalars_);
  ftmTreePP.setTreeType(ftm::TreeType::Join_Split);
  ftmTreePP.setVertexSoSoffsets(voffsets.data());
  ftmTreePP.setThreadNumber(threadNumber_);
  ftmTreePP.build<dataType, SimplexId>();
  ftmTreePP.setSegmentation(false);
  ftmTreePP.computePersistencePairs<dataType>(JTPairs, true);
  ftmTreePP.computePersistencePairs<dataType>(STPairs, false);

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::compressForPersistenceDiagram(
  int vertexNumber,
  dataType *inputData,
  dataType *outputData,
  const double &tol) {
  ttk::Timer t;
  ttk::Timer t1;

  std::vector<SimplexId> inputOffsets(vertexNumber);
  for(int i = 0; i < vertexNumber; ++i)
    inputOffsets[i] = i;

  // 1. Compute persistence pairs.
  std::vector<std::tuple<dataType, int>> topoIndices;

  // Compute global min & max
  dataType iter;
  dataType maxValue = inputData[0];
  int maxIndex = 0;
  for(int i = 1; i < vertexNumber; ++i) {
    iter = inputData[i];
    if(iter > maxValue) {
      maxIndex = i;
      maxValue = iter;
    }
  }

  dataType minValue = inputData[0];
  int minIndex = 0;
  for(int i = 1; i < vertexNumber; ++i) {
    iter = inputData[i];
    if(iter < minValue) {
      minIndex = i;
      minValue = iter;
    }
  }

  topoIndices.push_back(std::make_tuple(maxValue, maxIndex));
  topoIndices.push_back(std::make_tuple(minValue, minIndex));
  double tolerance = 0.01 * tol * (maxValue - minValue);
  double maxError = 0.01 * maximumError_ * (maxValue - minValue);

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Computed min/max in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  bool sqDomain = false;
  bool sqRange = false;

  const char *sq = sqMethod_.c_str();
  int nbCrit = 0;
  std::vector<int> simplifiedConstraints;
  if(strcmp(sq, "") == 0 && !zfpOnly_) {
    // No SQ: perform topological control

    std::vector<std::tuple<SimplexId, SimplexId, dataType>> JTPairs;
    std::vector<std::tuple<SimplexId, SimplexId, dataType>> STPairs;
    computePersistencePairs<dataType>(
      JTPairs, STPairs, inputData, inputOffsets.data());

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Persistence pairs computed in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }

    int nbJ = JTPairs.size();
    int nbS = STPairs.size();
    std::vector<int> critConstraints(2 * nbJ + 2 * nbS);

    topologicalSimplification.setupTriangulation(triangulation_);
    // auto* authorizedSaddles = new std::vector<int>();

    dataType maxEpsilon = 0;
    dataType epsilonSum2 = 0;
    dataType epsilonSum1 = 0;
    dataType persistentSum2 = 0;
    dataType persistentSum1 = 0;

    // Join
    for(int i = 0; i < nbJ; ++i) {
      SimplexId cp1 = std::get<0>(JTPairs[i]);
      SimplexId cp2 = std::get<1>(JTPairs[i]);
      dataType idt1 = inputData[cp1];
      dataType idt2 = inputData[cp2];
      dataType p1 = std::max(idt2, idt1) - std::min(idt2, idt1);
      if(p1 > tolerance) {
        persistentSum2 += (p1 * p1);
        persistentSum1 += abs<dataType>(p1);
        int type1 = topologicalSimplification.getCriticalType(
          cp1, inputData, inputOffsets.data());
        int type2 = topologicalSimplification.getCriticalType(
          cp2, inputData, inputOffsets.data());
        if(type1 == 0) {
          // authorizedSaddles->push_back(cp1);
          topoIndices.push_back(std::make_tuple(idt1, cp1));
        }
        if(type2 == 0) {
          // authorizedSaddles->push_back(cp2);
          topoIndices.push_back(std::make_tuple(idt2, cp2));
        }

        nbCrit += 2;
        critConstraints[2 * i] = cp1;
        critConstraints[2 * i + 1] = cp2;
      } else {
        if(maxEpsilon < abs<dataType>(p1))
          maxEpsilon = abs<dataType>(p1);
        epsilonSum2 += (p1 * p1);
        epsilonSum1 += abs<dataType>(p1);
        critConstraints[2 * i] = -1;
        critConstraints[2 * i + 1] = -1;
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Join pairs post-processed in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }

    // Split
    for(int i = nbJ; i < nbJ + nbS; ++i) {
      int si = i - nbJ;
      SimplexId cp1 = std::get<0>(STPairs[si]);
      SimplexId cp2 = std::get<1>(STPairs[si]);
      dataType idt1 = inputData[cp1];
      dataType idt2 = inputData[cp2];
      dataType p1 = std::max(idt2, idt1) - std::min(idt2, idt1);
      if(p1 > tolerance) {
        persistentSum2 += (p1 * p1);
        persistentSum1 += abs<dataType>(p1);
        // Saddle selection.
        int type1 = topologicalSimplification.getCriticalType(
          cp1, inputData, inputOffsets.data());
        int type2 = topologicalSimplification.getCriticalType(
          cp2, inputData, inputOffsets.data());
        if(type1 == 0) {
          // authorizedSaddles->push_back(cp1);
          topoIndices.push_back(std::make_tuple(idt1, cp1));
        }
        if(type2 == 0) {
          // authorizedSaddles->push_back(cp2);
          topoIndices.push_back(std::make_tuple(idt2, cp2));
        }

        nbCrit += 2;
        critConstraints[2 * i] = cp1;
        critConstraints[2 * i + 1] = cp2;
      } else {
        if(maxEpsilon < p1)
          maxEpsilon = abs<dataType>(p1);
        epsilonSum2 += (p1 * p1);
        epsilonSum1 += abs<dataType>(p1);
        critConstraints[2 * i] = -1;
        critConstraints[2 * i + 1] = -1;
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Split pairs post-processed in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }

    simplifiedConstraints.resize(nbCrit);
    {
      int j = 0;
      for(int i = 0; i < 2 * nbJ + 2 * nbS; ++i) {
        int c = critConstraints[i];
        if(c != -1)
          simplifiedConstraints[j++] = c;
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Got pairs in " << t.getElapsedTime()
          << " s. (" << threadNumber_ << " thread(s))." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }

    // 2. Perform topological simplification with constraints.
    if(useTopologicalSimplification_) {
      topologicalSimplification.setInputScalarFieldPointer(inputData);
      topologicalSimplification.setOutputScalarFieldPointer(outputData);
      topologicalSimplification.setInputOffsetScalarFieldPointer(
        inputOffsets.data());
      compressedOffsets_.resize(vertexNumber);
      for(int i = 0; i < vertexNumber; ++i)
        compressedOffsets_[i] = i;
      topologicalSimplification.setOutputOffsetScalarFieldPointer(
        compressedOffsets_.data());
      topologicalSimplification.setVertexIdentifierScalarFieldPointer(
        simplifiedConstraints.data());
      topologicalSimplification.setConstraintNumber(nbCrit);
      int status = 0;
      status = topologicalSimplification.execute<dataType, SimplexId>();
      if(status != 0) {
        return status;
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Performed simplification in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }

  } else if(strcmp(sq, "r") == 0 || strcmp(sq, "R") == 0) {
    // Range-based SQ: simple range quantization (automatic)
    sqRange = true;
  } else if(strcmp(sq, "d") == 0 || strcmp(sq, "D") == 0) {
    // Domain-based SQ: range quantization + later domain control
    sqDomain = true;
  } else if(zfpOnly_) {
    return 0;
  } else {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Unrecognized SQ option (" << sqMethod_
          << ")." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }

    return -3;
  }

  // 3. Subdivision and attribution of a topological index

  auto cmp = [](const std::tuple<dataType, int> &a,
                const std::tuple<dataType, int> &b) {
    return std::get<0>(a) < std::get<0>(b);
  };
  std::sort(topoIndices.begin(), topoIndices.end(), cmp);
  std::vector<std::tuple<dataType, int>> segments;
  // ith rank = ith segment (from bottom)
  // 0: value
  // 1: critical point index or -1 if !critical
  bool subdivide = !dontSubdivide_;
  auto l = (int)topoIndices.size();
  if(l < 1) {
    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Trivial subdivision performed."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
    segments.push_back(topoIndices[l - 1]);
  } else {
    for(int i = 0; i < l; ++i) {
      dataType v0 = std::get<0>(topoIndices[i]);
      if(i == l - 1) {
        segments.push_back(topoIndices[i]);
        break;
      }

      dataType v1 = std::get<0>(topoIndices[i + 1]);
      int i1 = std::get<1>(topoIndices[i + 1]);
      auto diff = (double)(v1 - v0);

      if(diff == 0)
        continue;
      if(!subdivide || (diff < (dataType)maxError)) {
        segments.push_back(std::make_tuple(v1, i1));
      } else {
        // Subdivide.
        double nSegments = std::ceil(diff / maxError);
        for(int j = 0, nbs = (int)nSegments; j < nbs; ++j) {
          dataType sample = v0 + j * maxError;
          int int1 = i1;
          segments.push_back(std::make_tuple(sample, int1));
        }
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Subdivision/topological index"
        << " attribution in " << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  // 4. Affect segment value to all points.
  segmentation_.resize(vertexNumber);
  std::sort(segments.begin(), segments.end(), cmp);
  for(int i = 0; i < vertexNumber; ++i) {
    dataType scalar = inputData[i];

    auto begin = segments.begin();
    auto end = segments.end();
    auto it = std::lower_bound(
      segments.begin(), segments.end(), std::make_tuple(scalar, -1), cmp);

    if(it != end) {
      std::tuple<dataType, int> tt = *it;
      int j = it - begin;
      int last = (int)segments.size() - 1;
      if(j < last) {
        dataType dtv = std::get<0>(tt);
        if(j > 0) {
          dtv = 0.5 * (dtv + std::get<0>(*(it - 1)));
        }

        outputData[i] = dtv;
        int seg = j;
        segmentation_[i] = seg;
      } else {
        segmentation_[i] = last;
        outputData[i] = std::get<0>(tt); // maxValue;
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Affected values in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  // 4.1. Simplify mapping
  auto segmentsSize = (int)segments.size();
  std::vector<bool> affectedSegments;
  affectedSegments.resize(segmentsSize);
  std::vector<int> oob;
  for(int i = 0; i < vertexNumber; ++i) {
    int seg = segmentation_[i];
    if(seg >= segmentsSize
       && std::find(oob.begin(), oob.end(), seg) == oob.end())
      oob.push_back(seg);
    else if(seg >= 0)
      affectedSegments[seg] = true;
    else {
      std::stringstream msg;
      msg << "[TopologicalCompression] Negative segment encoutered (" << i
          << ")" << std::endl;
      dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }
  }
  std::vector<int> empty;
  for(int i = 0; i < segmentsSize; ++i) {
    if(!affectedSegments[i])
      empty.push_back(i);
  }

  std::sort(oob.begin(), oob.end());
  std::sort(empty.begin(), empty.end());
  int indexLast = -1;
  if(oob.size() > empty.size()) {
    std::stringstream msg;
    msg << "[TopologicalCompression] WARN oob size > empty size: " << oob.size()
        << ", " << empty.size() << std::endl;
    dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  } else {
    // Replace
    for(int i = 0; i < vertexNumber; ++i) {
      int seg = segmentation_[i];
      if(seg >= segmentsSize) {
        auto begin = oob.begin();
        auto end = oob.end();
        auto it = std::lower_bound(begin, end, seg);
        if(it != end) {
          int j = (int)(it - begin);
          segmentation_[i] = empty[j];
          affectedSegments[j] = true;
        }
      } else if(!affectedSegments[seg]) {
        std::stringstream msg;
        msg << "[TopologicalCompression] Something impossible happened"
            << std::endl;
        dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
      }
    }

    if(empty.size() > oob.size()) {
      std::vector<int> map2(segmentsSize, -1);
      for(int i = 0; i < segmentsSize; ++i) {
        if(affectedSegments[i])
          continue;
        for(int j = (segmentsSize - 1); j > i; --j) {
          if(!affectedSegments[j] || map2[j] > 0)
            continue;
          map2[j] = i;
          affectedSegments[j] = false;
          affectedSegments[i] = true;
          break;
        }
      }

      bool doneAffecting = false;
      for(int i = 0; i < segmentsSize; ++i) {
        if(!affectedSegments[i]) {
          doneAffecting = true;
        } else {
          if(doneAffecting) {
            std::stringstream msg;
            msg << "[TopologicalCompression] Hole detected at " << i
                << "th segment." << std::endl;
            dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
          } else {
            indexLast = i;
          }
        }
      }

      for(int i = 0; i < vertexNumber; ++i) {
        int seg = segmentation_[i];
        if(map2[seg] > 0)
          segmentation_[i] = map2[seg];
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Simplified mapping in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  // 5. Expose mapping.
  std::vector<bool> already(vertexNumber);
  for(int i = 0; i < vertexNumber; ++i)
    already[i] = false; // init

  for(int i = 0; i < vertexNumber; ++i) {
    int vert = segmentation_[i];
    if(!already[vert]) {
      already[vert] = true;
      dataType dttt = outputData[i];
      mapping_.push_back(std::make_tuple((double)dttt, vert));
    }
  }

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Exposed mapping in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  // 6. SQ-D correction step.
  // Indices stay compact.
  if(sqDomain && !sqRange) {
    std::vector<bool> markedVertices(vertexNumber); // false by default
    std::vector<bool> markedSegments(segmentsSize); // false by default
    int newIndex = segmentsSize;

    for(int i = 0; i < vertexNumber; ++i) {
      // Skip processed vertices.
      if(markedVertices[i])
        continue;

      int seg = segmentation_[i];
      bool newSegment = markedSegments[seg];
      dataType minNewSegment = inputData[i];
      dataType maxNewSegment = minNewSegment;

      if(!newSegment) {
        markedSegments[seg] = true;
      }

      // Neighborhood search with seed = i
      std::stack<int> s;
      s.push(i);

      while(!s.empty()) {
        // Get next element.
        int vertex = s.top();
        s.pop();

        // Mark vertex as processed.
        markedVertices[vertex] = true;

        // New region was detected.
        if(newSegment) {
          // Affect new segmentation id to all vertices in region.
          segmentation_[vertex] = newIndex;
          // Progressively compute min & max on local region.
          dataType currentValue = inputData[vertex];
          if(currentValue > maxNewSegment)
            maxNewSegment = currentValue;
          if(currentValue < minNewSegment)
            minNewSegment = currentValue;
        }

        // Get neighbors.
        SimplexId neighborNumber
          = triangulation_->getVertexNeighborNumber(vertex);
        for(SimplexId j = 0; j < neighborNumber; ++j) {
          SimplexId neighbor;
          triangulation_->getVertexNeighbor(vertex, j, neighbor);

          // Add current neighbor to processing stack.
          if(!markedVertices[neighbor] && segmentation_[neighbor] == seg) {
            s.push(neighbor);
          }
        }
      }

      // Affect a value to the new segmentation index.
      if(newSegment) {
        dataType result = (maxNewSegment + minNewSegment) / 2;
        segments.push_back(std::make_tuple(result, newIndex));
        mapping_.push_back(std::make_tuple((double)result, newIndex));
        newIndex++; // Set index for next potential segment.
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] SQ-D correction step in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }
  }

  // 7. [ZFP]: max constraints, min constraints
  if(!sqDomain && !sqRange && !zfpOnly_) {
    for(int i = 0; i < nbCrit; ++i) {
      SimplexId id = simplifiedConstraints[i];
      dataType val = inputData[id];
      int type = topologicalSimplification.getCriticalType(
        id, inputData, inputOffsets.data());
      if(type == -1 // Local_minimum
         || type == 1 // Local_maximum
         || type == 0) {
        criticalConstraints_.push_back(std::make_tuple(id, (double)val, type));
      } else { // 0 -> Saddle
      }
    }

    {
      std::stringstream msg;
      msg << "[TopologicalCompression] Exposed constraints in "
          << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      t.reStart();
    }
  }

  {
    if(indexLast > -1) {
      if(indexLast + 1 != (int)mapping_.size()) {
        std::stringstream msg;
        msg << "[TopologicalCompression] possible affectation mismatch "
            << "(" << (indexLast + 1) << ", " << mapping_.size() << ")"
            << std::endl;
        dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
      }
    }

    int nSegments = indexLast > -1 ? indexLast + 1 : (int)segments.size() - 1;
    this->nbSegments = nSegments;
    this->nbVertices = vertexNumber;
    std::stringstream msg;
    msg << "[TopologicalCompression] Affected " << nSegments << " segment"
        << (nbSegments > 1 ? "s" : "") << "." << std::endl;
    msg << "[TopologicalCompression] Data-set (" << vertexNumber
        << " points) processed in " << t1.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // TTK_PERSISTENCEDIAGRAMCOMPRESSION_H
