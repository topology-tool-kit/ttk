/// \ingroup base
/// \class ttk::BottleneckDistance 
/// \author Maxime Soler <soler.maxime@total.com>
/// \date The Date Here.
///
/// \brief TTK %bottleneckDistance processing package.
///
/// %BottleneckDistance is a TTK processing package that 
/// takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkBottleneckDistance.cpp %for a usage example.

#ifndef _BOTTLENECKDISTANCE_H
#define _BOTTLENECKDISTANCE_H

#ifndef diagramTuple
#define diagramTuple tuple<idVertex, NodeType, idVertex, NodeType, dataType, idVertex, \
  dataType, float, float, float, dataType, float, float, float>
#endif

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include                  <Munkres.h>

#include                  <string>

namespace ttk {

  class BottleneckDistance : public Debug {

    public:
        
      BottleneckDistance();
      
      ~BottleneckDistance();

      /// Execute the package.
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh 
      /// traversals, the function setupTriangulation() must be called on this 
      /// object prior to this function, in a clearly distinct pre-processing 
      /// steps. An error will be returned otherwise.
      /// \note In such a case, it is recommended to exclude 
      /// setupTriangulation() from any time performance measurement.
      /// \param argment Dummy integer argument.
      /// \return Returns 0 upon success, negative values otherwise.
      template <typename dataType>
        int execute(bool usePersistenceMetric, double alpha);
    
      template <typename dataType>
        int computeBottleneck(
        vector<diagramTuple> *CTDiagram1,
        vector<diagramTuple> *CTDiagram2,
        vector<tuple<idVertex, idVertex, dataType> > *matchings,
        bool usePersistenceMetric,
        double alpha);

      inline int setCTDiagram1(void *diagram) {
        outputCT1_ = diagram;
        return 0;
      }

      inline int setCTDiagram2(void *diagram) {
        outputCT2_ = diagram;
        return 0;
      }

      inline int setOutputMatchings(void* matchings) {
        matchings_ = matchings;
        return 0;
      }

      inline int setWasserstein(const string &wasserstein) {
        wasserstein_ = wasserstein;
        return 0;
      }

      template <typename dataType>
        dataType
        getDistance();

      template<typename type>
      static type abs(const type var){
         return (var >= 0) ? var : -var;
      }

      template<typename type>
      static type abs_diff(const type var1, const type var2){
         return (var1 > var2) ? var1 - var2 : var2 - var1;
      }

      // Triangulation
      inline int setupTriangulation(Triangulation *triangulation) 
      {
        triangulation_ = triangulation;
        if (triangulation_) triangulation_->preprocessVertexNeighbors();        
        return 0;
      }
    
    protected:
    
      void                      *outputCT1_;
      void                      *outputCT2_;
      void                      *matchings_; // ids from CT1 to CT2
      void                      *distance_;

      string                    wasserstein_;
      Triangulation             *triangulation_;
    
  private:

    template <typename dataType>
      bool isValidMatching(
        const vector<tuple<idVertex, idVertex, dataType>>* matchings,
        dataType thresholdMin) const;
  };
}

template <typename dataType>
  dataType
  BottleneckDistance::getDistance()
{
    return *static_cast<dataType*> (distance_);
}

template <typename dataType>
  int BottleneckDistance::execute(bool usePersistenceMetric, double alpha)
{
  
  Timer t;
  
  computeBottleneck(
    static_cast<vector<diagramTuple>*> (outputCT1_),
    static_cast<vector<diagramTuple>*> (outputCT2_),
    static_cast<vector<tuple<idVertex, idVertex, dataType>>*> (matchings_),
    usePersistenceMetric,
    alpha);

  {
    stringstream msg;
    msg << "[BottleneckDistance] Data-set processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    msg << "[BottleneckDistance] distance = " << *(dataType*) distance_ << endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <typename dataType>
int BottleneckDistance::computeBottleneck(
//  dataType* inputData1,
//  dataType* inputData2,
  // vector <   -- diagram
  // tuple <    -- pair of critical points
  //  idVertex
  //  NodeType
  //  idVertex
  //  NodeType
  //  dataType  -- persistance of pair
  //  idVertex  -- type (0/min, 1/saddle, 2/max)
    // dataType            -- scalar value at vertex 1
    // float, float, float -- vertex 1 coordinates
  vector<diagramTuple> *CTDiagram1,
  vector<diagramTuple> *CTDiagram2,
  vector<tuple<idVertex, idVertex, dataType>>* matchings,
  bool usePersistenceMetric,
  double alpha)
{
  dataType* distance = new dataType;

  int d1Size = CTDiagram1->size();
  int d2Size = CTDiagram2->size();

  // Discard invalid alpha parameter.
  if (alpha < 0.0 || alpha > 1.0)
    alpha = 1.0;

  // Compute geometrical range.
  float minX1, maxX1, minY1, maxY1, minZ1, maxZ1;
  float minX2, maxX2, minY2, maxY2, minZ2, maxZ2;
  float minX, minY, minZ, maxX, maxY, maxZ;
  minX1 = minY1 = minZ1 = minX2 = minY2 = minZ2 = std::numeric_limits<float>::max();
  maxX1 = maxY1 = maxZ1 = maxX2 = maxY2 = maxZ2 = std::numeric_limits<float>::min();

  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    float xa = get<7>(t), ya = get<8>(t), za = get<9>(t);
    float xb = get<11>(t), yb = get<12>(t), zb = get<13>(t);
    minX1 = std::min(std::min(minX1, xa), xb);
    minY1 = std::min(std::min(minY1, ya), yb);
    minZ1 = std::min(std::min(minZ1, za), zb);
    maxX1 = std::max(std::max(maxX1, xa), xb);
    maxY1 = std::max(std::max(maxY1, ya), yb);
    maxZ1 = std::max(std::max(maxZ1, za), zb);
  }

  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    float xa = get<7>(t), ya = get<8>(t), za = get<9>(t);
    float xb = get<11>(t), yb = get<12>(t), zb = get<13>(t);
    minX2 = std::min(std::min(minX2, xa), xb);
    minY2 = std::min(std::min(minY2, ya), yb);
    minZ2 = std::min(std::min(minZ2, za), zb);
    maxX2 = std::max(std::max(maxX2, xa), xb);
    maxY2 = std::max(std::max(maxY2, ya), yb);
    maxZ2 = std::max(std::max(maxZ2, za), zb);
  }

  minX = std::min(minX1, minX2); maxX = std::max(maxX1, maxX2);
  minY = std::min(minY1, minY2); maxY = std::max(maxY1, maxY2);
  minZ = std::min(minZ1, minZ2); maxZ = std::max(maxZ1, maxZ2);

  double maxDistance = std::sqrt(
      std::pow(maxX - minX, 2) + std::pow(maxY - minY, 2) + std::pow(maxZ - minZ, 2));

  // Limit computation time.
  vector<dataType> toSort;
  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    dataType persistence = abs<dataType>(get<4>(t));
    toSort.push_back(persistence);
  }
  for (int i = 0; i < d2Size; ++i) {
    diagramTuple t = CTDiagram2->at(i);
    dataType persistence = abs<dataType>(get<4>(t));
    toSort.push_back(persistence);
  }
  std::sort(toSort.begin(), toSort.end());
  double epsilon = 0.0000001;
  int largeSize = 2000;
  dataType zeroThresh = (dataType) epsilon;
  if (d1Size + d2Size > largeSize + 1) {
    zeroThresh = toSort[d1Size + d2Size - largeSize];
    if (toSort[d1Size + d2Size - (largeSize+1)] == zeroThresh)
      zeroThresh += (dataType) epsilon;
  }
  if (zeroThresh < epsilon) zeroThresh = epsilon;

  vector<tuple<idVertex, idVertex, dataType>>* minMatchings
    = new vector<tuple<idVertex, idVertex, dataType>>();
  vector<tuple<idVertex, idVertex, dataType>>* maxMatchings
    = new vector<tuple<idVertex, idVertex, dataType>>();
  vector<tuple<idVertex, idVertex, dataType>>* sadMatchings
      = new vector<tuple<idVertex, idVertex, dataType>>();

  Munkres* solverMin = new Munkres();
  Munkres* solverMax = new Munkres();
  Munkres* solverSad = new Munkres();

  // Init cost matrices.
  int wasserstein = -1;
  if (wasserstein_ != "inf") {
    int n = stoi(wasserstein_);
    if (n < 1) return -4;
    wasserstein = n;
  }

  int nbRowMin = 0, nbColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int nbRowSad = 0, nbColSad = 0;

  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    NodeType nt1 = get<1>(t);
    NodeType nt2 = get<3>(t);
    dataType dt = get<4>(t);
    if (abs<dataType>(dt) < zeroThresh) continue;

    if (nt1 == NodeType::Local_minimum && nt2 == NodeType::Local_maximum) ++nbRowMax;
    else {
      if (nt1 == NodeType::Local_maximum || nt2 == NodeType::Local_maximum) ++nbRowMax;
      if (nt1 == NodeType::Local_minimum || nt2 == NodeType::Local_minimum) ++nbRowMin;
      if ((nt1 == NodeType::Saddle1 && nt2 == NodeType::Saddle2)
        || (nt1 == NodeType::Saddle2 && nt2 == NodeType::Saddle1)) ++nbRowSad;
    }
  }

  for (int i = 0; i < d2Size; ++i) {
    diagramTuple t = CTDiagram2->at(i);
    NodeType nt1 = get<1>(t);
    NodeType nt2 = get<3>(t);
    dataType dt = get<4>(t);
    if (abs<dataType>(dt) < zeroThresh) continue;

    if (nt1 == NodeType::Local_minimum && nt2 == NodeType::Local_maximum) ++nbColMax;
    else {
      if (nt1 == NodeType::Local_maximum || nt2 == NodeType::Local_maximum) ++nbColMax;
      if (nt1 == NodeType::Local_minimum || nt2 == NodeType::Local_minimum) ++nbColMin;
      if ((nt1 == NodeType::Saddle1 && nt2 == NodeType::Saddle2)
        || (nt1 == NodeType::Saddle2 && nt2 == NodeType::Saddle1)) ++nbColSad;
    }
  }

  // Fill cost matrices.
  int nMin = std::max(nbColMin, nbRowMin);
  int nMax = std::max(nbColMax, nbRowMax);
  int nSad = std::max(nbColSad, nbRowSad);

  dataType** minMatrix = new dataType*[nMin];
  dataType** maxMatrix = new dataType*[nMax];
  dataType** sadMatrix = new dataType*[nSad];

  for (int i = 0; i < nMin; ++i) minMatrix[i] = new dataType[nMin];
  for (int i = 0; i < nMax; ++i) maxMatrix[i] = new dataType[nMax];
  for (int i = 0; i < nSad; ++i) sadMatrix[i] = new dataType[nSad];

  int maxI = 0, minI = 0,
      maxJ = 0, minJ = 0,
      sadI = 0, sadJ = 0;

  for (int i = 0; i < d1Size; ++i)
  {
    diagramTuple t1 = CTDiagram1->at(i);
    if (abs<dataType>(get<4>(t1)) < zeroThresh) continue;

    bool isMin1 = (get<1>(t1) == NodeType::Local_minimum || get<3>(t1) == NodeType::Local_minimum);
    bool isMax1 = (get<1>(t1) == NodeType::Local_maximum || get<3>(t1) == NodeType::Local_maximum);
    bool isSad1 = ((get<1>(t1) == NodeType::Saddle1 && get<3>(t1) == NodeType::Saddle2) ||
        (get<1>(t1) == NodeType::Saddle2 && get<3>(t1) == NodeType::Saddle1));

    if (get<1>(t1) == NodeType::Local_minimum && get<3>(t1) == NodeType::Local_maximum) {
      isMin1 = false;
      isMax1 = true;
    }
    dataType rX = get<6>(t1);
    dataType rY = get<10>(t1);
    dataType rDiff2 = usePersistenceMetric ?
                      wasserstein > 0 ?
                        std::pow(abs<dataType>(get<4>(t1)), wasserstein) :
                        abs<dataType>(get<4>(t1)) : 0;

    // Reinit indices.
    minJ = 0;
    maxJ = 0;
    sadJ = 0;

    for (int j = 0; j < d2Size; ++j)
    {
      diagramTuple t2 = CTDiagram2->at(j);
      if (abs<dataType>(get<4>(t2)) < zeroThresh) continue;

      bool isMin2 = (get<1>(t2) == NodeType::Local_minimum || get<3>(t2) == NodeType::Local_minimum);
      bool isMax2 = (get<1>(t2) == NodeType::Local_maximum || get<3>(t2) == NodeType::Local_maximum);
      bool isSad2 = ((get<1>(t2) == NodeType::Saddle1 && get<3>(t2) == NodeType::Saddle2) ||
                     (get<1>(t2) == NodeType::Saddle2 && get<3>(t2) == NodeType::Saddle1));

      if (get<1>(t2) == NodeType::Local_minimum && get<3>(t2) == NodeType::Local_maximum) {
        isMin2 = false;
        isMax2 = true;
      }
      if ((isMin1 && !isMin2) || (isMax1 && !isMax2) || (isSad1 && !isSad2)) continue;

      dataType cX = get<6>(t2);
      dataType cY = get<10>(t2);
      dataType cDiff2 = usePersistenceMetric ?
        wasserstein > 0 ?
          std::pow(abs<dataType>(get<4>(t2)), wasserstein) :
          abs<dataType>(get<4>(t2)) : 0;

      // 1. Compute Wasserstein distance (w = 2 => square mean; w = infty => bottleneck distance).
      // 2. Account for non-paired cost by removing the feature heights.
      dataType x = abs_diff<dataType>(rX, cX);
      dataType y = abs_diff<dataType>(rY, cY);

      double dist = std::sqrt(
          std::pow((get<7>(t1)+get<11>(t1))/2 - (get<7>(t2)+get<11>(t2))/2, 2) +
          std::pow((get<8>(t1)+get<12>(t1))/2 - (get<8>(t2)+get<12>(t2))/2, 2) +
          std::pow((get<9>(t1)+get<13>(t1))/2 - (get<9>(t2)+get<13>(t2))/2, 2)
      );
      dist /= maxDistance;

      if (isMin1 && isMin2) {
        double val = (wasserstein > 0 ? // Wasserstein
          std::pow(x, wasserstein) + std::pow(y, wasserstein) : // Bottleneck
          std::max(x, y)) - rDiff2 - cDiff2;
        minMatrix[minI][minJ++] = alpha * val + (1.0 - alpha) * dist;
      }
      else if (isMax1 && isMax2) {
        double val = (wasserstein > 0 ?
          std::pow(x, wasserstein) + std::pow(y, wasserstein) :
          std::max(x, y)) - rDiff2 - cDiff2;
        maxMatrix[maxI][maxJ++] = alpha * val + (1.0 - alpha) * dist;
      }
      else if (isSad1 && isSad2) {
        double val = (wasserstein > 0 ?
          std::pow(x, wasserstein) + std::pow(y, wasserstein) :
          std::max(x, y)) - rDiff2 - cDiff2;
        sadMatrix[sadI][sadJ++] = alpha * val + (1.0 - alpha) * dist;
      }

      // matrix[i][j] = x*x + y*y - rDiff2 - cDiff2; // (cY-cX)²
    }

    if (isMin1) ++minI;
    if (isMax1) ++maxI;
    if (isSad1) ++sadI;
  }


  if (d1Size > d2Size) {
    minI = 0;
    maxI = 0;
    sadI = 0;
    for (int i = 0; i < d1Size; ++i) {
      diagramTuple t1 = CTDiagram1->at(i);
      if (abs<dataType>(get<4>(t1)) < zeroThresh) continue;
      bool isMin1 = (get<1>(t1) == NodeType::Local_minimum || get<3>(t1) == NodeType::Local_minimum);
      bool isMax1 = (get<1>(t1) == NodeType::Local_maximum || get<3>(t1) == NodeType::Local_maximum);
      bool isSad1 = ((get<1>(t1) == NodeType::Saddle1 && get<3>(t1) == NodeType::Saddle2) ||
                     (get<1>(t1) == NodeType::Saddle2 && get<3>(t1) == NodeType::Saddle1));

      if (get<1>(t1) == NodeType::Local_minimum && get<3>(t1) == NodeType::Local_maximum) {
        isMin1 = false;
        isMax1 = true;
      }
      if (!isMin1 && !isMax1 && !isSad1) continue;

      dataType rDiff2 = wasserstein > 0 ?
                        std::pow(abs<dataType>(get<4>(t1)), wasserstein) :
                        abs<dataType>(get<4>(t1));

      if (isMin1) {
        int start = std::min(nbColMin, nbRowMin);
        int end = std::max(nbColMin, nbRowMin);

        if (nbColMin < nbRowMin) {
          for (int j = start; j < end; ++j)
            minMatrix[minI][j] = rDiff2;
        } else {
          for (int j = start; j < end; ++j)
            minMatrix[j][minI] = rDiff2;
        }

        ++minI;
      }
      else if (isMax1) {
        int start = std::min(nbColMax, nbRowMax);
        int end = std::max(nbColMax, nbRowMax);

        if (nbColMax < nbRowMax) {
          for (int j = start; j < end; ++j)
            maxMatrix[maxI][j] = rDiff2;
        } else {
          for (int j = start; j < end; ++j)
            maxMatrix[j][maxI] = rDiff2;
        }

        ++maxI;
      }
      else if (isSad1) {
        int start = std::min(nbColSad, nbRowSad);
        int end = std::max(nbColSad, nbRowSad);

        if (nbColSad < nbRowSad) {
          for (int j = start; j < end; ++j)
            sadMatrix[sadI][j] = rDiff2;
        } else {
          for (int j = start; j < end; ++j)
            sadMatrix[j][sadI] = rDiff2;
        }

        ++sadI;
      }
    }
  } else if (d2Size > d1Size) {
    minJ = 0;
    maxJ = 0;
    for (int i = 0; i < d2Size; ++i) {
      diagramTuple t2 = CTDiagram2->at(i);
      if (abs<dataType>(get<4>(t2)) < zeroThresh) continue;
      bool isMin2 = (get<1>(t2) == NodeType::Local_minimum || get<3>(t2) == NodeType::Local_minimum);
      bool isMax2 = (get<1>(t2) == NodeType::Local_maximum || get<3>(t2) == NodeType::Local_maximum);
      bool isSad2 = ((get<1>(t2) == NodeType::Saddle1 && get<3>(t2) == NodeType::Saddle2) ||
                     (get<1>(t2) == NodeType::Saddle2 && get<3>(t2) == NodeType::Saddle1));

      if (get<1>(t2) == NodeType::Local_minimum && get<3>(t2) == NodeType::Local_maximum) {
        isMin2 = false;
        isMax2 = true;
      }

      if (!isMin2 && !isMax2 && !isSad2) continue;

      dataType cDiff2 = wasserstein > 0 ?
                        std::pow(abs<dataType>(get<4>(t2)), wasserstein) :
                        abs<dataType>(get<4>(t2));

      if (isMin2) {
        int start = std::min(nbColMin, nbRowMin);
        int end = std::max(nbColMin, nbRowMin);

        if (nbColMin < nbRowMin) {
          for (int j = start; j < end; ++j)
            minMatrix[minJ][j] = cDiff2;
        } else {
          for (int j = start; j < end; ++j)
            minMatrix[j][minJ] = cDiff2;
        }

        ++minJ;
      }
      else if (isMax2) {
        int start = std::min(nbColMax, nbRowMax);
        int end = std::max(nbColMax, nbRowMax);

        if (nbColMax < nbRowMax) {
          for (int j = start; j < end; ++j)
            maxMatrix[maxJ][j] = cDiff2;
        } else {
          for (int j = start; j < end; ++j)
            maxMatrix[j][maxJ] = cDiff2;
        }

        ++maxJ;
      }
      else if (isSad2) {
        int start = std::min(nbColSad, nbRowSad);
        int end = std::max(nbColSad, nbRowSad);

        if (nbColSad < nbRowSad) {
          for (int j = start; j < end; ++j)
            sadMatrix[sadJ][j] = cDiff2;
        } else {
          for (int j = start; j < end; ++j)
            sadMatrix[j][sadJ] = cDiff2;
        }

        ++sadJ;
      }
    }
  }

  if (wasserstein > 0) {
    solverMin->setInput(nMin, nMin, (void*) minMatrix);
    solverMin->run<dataType>(minMatchings);

    solverMax->setInput(nMax, nMax, (void*) maxMatrix);
    solverMax->run<dataType>(maxMatchings);

    solverSad->setInput(nSad, nSad, (void*) sadMatrix);
    solverSad->run<dataType>(sadMatchings);
  } else {

    // Cf. "A New Algorithm for Solving Linear Bottleneck Assignment Problem"
    // P. S. Pundir, S. K. Porwal, B. P. Singh
    // Journal of Institute of Science and Technology.

    // Preprocess bottleneck.
    // 1. for minima
    dataType minCol, minRow;
    vector<dataType> minRowCols;
    minRowCols.reserve(nbRowMin + nbColMin);
    for (int i = 0; i < nbRowMin; ++i) {
      minRow = minMatrix[i][0];
      for (int j = 1; j < nbColMin; ++j) {
        minRow = std::min(minRow, minMatrix[i][j]);
      }
      minRowCols[i] = minRow;
    }
    for (int j = 0; j < nbColMin; ++j) {
      minCol = minMatrix[0][j];
      for (int i = 0; i < nbRowMin; ++i) {
        minCol = std::min(minCol, minMatrix[i][j]);
      }
      minRowCols[nbRowMin + j] = minCol;
    }
    dataType thresholdMin = minRowCols[0];
    for (int i = 1, s = minRowCols.size(); i < s; ++s)
      thresholdMin = std::max(thresholdMin, minRowCols[i]);

    // 2. for maxima
    dataType maxCol, maxRow;
    vector<dataType> maxRowCols;
    maxRowCols.reserve(nbRowMax + nbColMax);
    dataType thresholdMax = maxMatrix[0][0];
    for (int i = 0; i < nbRowMax; ++i) {
      maxRow = maxMatrix[i][0];
      for (int j = 1; j < nbColMax; ++j) {
        maxRow = std::min(maxRow, maxMatrix[i][j]);
      } maxRowCols.push_back(maxRow);
    }
    for (int j = 1; j < nbColMax; ++j) {
      maxCol = maxMatrix[0][j];
      for (int i = 0; i < nbRowMax; ++i) {
        maxCol = std::min(maxCol, maxMatrix[i][j]);
      } maxRowCols.push_back(maxCol);
    }
    for (int i = 0, s = maxRowCols.size(); i < s; ++i)
      thresholdMax = std::max(thresholdMax, maxRowCols[i]);

    // 3. for saddles
    dataType sadCol, sadRow;
    vector<dataType> sadRowCols;
    sadRowCols.reserve(nbRowSad + nbColSad);
    dataType thresholdSad = sadMatrix[0][0];
    for (int i = 0; i < nbRowSad; ++i) {
      sadRow = sadMatrix[i][0];
      for (int j = 1; j < nbColSad; ++j) {
        sadRow = std::min(sadRow, sadMatrix[i][j]);
      } sadRowCols.push_back(sadRow);
    }
    for (int j = 1; j < nbColSad; ++j) {
      sadCol = sadMatrix[0][j];
      for (int i = 1; i < nbRowSad; ++i) {
        sadCol = std::min(sadCol, sadMatrix[i][j]);
      } sadRowCols.push_back(sadCol);
    }
    for (int i = 0, s = sadRowCols.size(); i < s; ++i)
      thresholdSad = std::max(thresholdSad, sadRowCols[i]);

    // Reorg.
    dataType** minBottleneckMatrix = new dataType*[nMin];
    for (int i = 0; i < nMin; ++i)  {
      minBottleneckMatrix[i] = new dataType[nMin];

      for (int j = 0; j < nMin; ++j) {
        if (i >= nbRowMin ||j >= nbColMin)
          minBottleneckMatrix[i][j] = minMatrix[i][j];
        if (i < nbRowMin && j < nbColMin) {
          minBottleneckMatrix[i][j] = minMatrix[i][j] > thresholdMin ?
            std::numeric_limits<dataType>::max() : minMatrix[i][j];
        }
      }
    }

    dataType** maxBottleneckMatrix = new dataType*[nMax];
    for (int i = 0; i < nMax; ++i) {
      maxBottleneckMatrix[i] = new dataType[nMax];

      for (int j = 0; j < nMax; ++j) {
        if (i >= nbRowMax || j >= nbColMax)
          maxBottleneckMatrix[i][j] = maxMatrix[i][j];
        if (i < nbRowMax && j < nbColMax) {
          maxBottleneckMatrix[i][j] = maxMatrix[i][j] > thresholdMax ?
            std::numeric_limits<dataType>::max() : maxMatrix[i][j];
        }
      }
    }

    dataType** sadBottleneckMatrix = new dataType*[nSad];
    for (int i = 0; i < nSad; ++i) {
      sadBottleneckMatrix[i] = new dataType[nSad];

      for (int j = 0; j < nSad; ++j) {
        if (i >= nbRowSad || j >= nbColSad)
          sadBottleneckMatrix[i][j] = sadMatrix[i][j];
        if (i < nbRowSad && j < nbColSad) {
          sadBottleneckMatrix[i][j] = sadMatrix[i][j] > thresholdSad ?
            std::numeric_limits<dataType>::max() : sadMatrix[i][j];
        }
      }
    }

    // Apply Munkres for the first time.

    solverMin->setInput(nMin, nMin, (void*) minBottleneckMatrix);
    solverMin->run<dataType>(minMatchings);

    solverMax->setInput(nMax, nMax, (void*) maxBottleneckMatrix);
    solverMax->run<dataType>(maxMatchings);

    solverSad->setInput(nSad, nSad, (void*) sadBottleneckMatrix);
    solverSad->run<dataType>(sadMatchings);

    // Check assignment maximum (infty => restart Munkres with a higher - next - threshold).
    int minIter = 2;
    while (!isValidMatching(minMatchings, thresholdMin) && minIter <= nbRowMin * nbColMin) {
      minMatchings->clear();
      // Find next threshold.
      dataType newThresholdMin = std::numeric_limits<dataType>::max();
      for (int i = 0; i < nbRowMin; ++i) {
        for (int j = 0; j < nbColMin; ++j) {
          dataType e = minMatrix[i][j];
          if (e > thresholdMin && e < newThresholdMin) newThresholdMin = e;
        }
      }
      thresholdMin = newThresholdMin;

      // Apply Munkres again.
      for (int i = 0; i < nbRowMin; ++i) {
        for (int j = 0; j < nbColMin; ++j) {
          dataType e = minMatrix[i][j];
          minBottleneckMatrix[i][j] = (e > thresholdMin) ?
            std::numeric_limits<dataType>::max() : e;
        }
      }

      solverMin->setInput(nMin, nMin, (void*) minBottleneckMatrix);
      solverMin->run<dataType>(minMatchings);
      minIter++;
    }

    if (minBottleneckMatrix != NULL) {
      for (int r = 0; r < nMin; ++r)
        delete[] (minBottleneckMatrix[r]);
      delete[] minBottleneckMatrix;
      minBottleneckMatrix = NULL;
    }

    int maxIter = 2;
    while (!isValidMatching(maxMatchings, thresholdMax) && maxIter <= nbRowMax * nbColMax) {
      maxMatchings->clear();
      // Find next threshold.
      dataType newThresholdMax = std::numeric_limits<dataType>::max();
      for (int i = 0; i < nbRowMax; ++i) {
        for (int j = 0; j < nbColMax; ++j) {
          dataType e = maxMatrix[i][j];
          if (e > thresholdMax && e < newThresholdMax) newThresholdMax = e;
        }
      }
      thresholdMax = newThresholdMax;

      // Apply Munkres again.
      for (int i = 0; i < nbRowMax; ++i) {
        for (int j = 0; j < nbColMax; ++j) {
          dataType e = maxMatrix[i][j];
          maxBottleneckMatrix[i][j] = (e > thresholdMax) ?
            std::numeric_limits<dataType>::max() : e;
        }
      }

      solverMax->setInput(nMax, nMax, (void*) maxBottleneckMatrix);
      solverMax->run<dataType>(maxMatchings);
      maxIter++;
    }

    if (maxBottleneckMatrix != NULL) {
      for (int r = 0; r < nMax; ++r)
        delete[] (maxBottleneckMatrix[r]);
      delete[] maxBottleneckMatrix;
      maxBottleneckMatrix = NULL;
    }

    int sadIter = 2;
    while (!isValidMatching(sadMatchings, thresholdSad) && sadIter <= nbRowSad * nbColSad) {
      sadMatchings->clear();
      // Find next threshold.
      dataType newThresholdSad = std::numeric_limits<dataType>::max();
      for (int i = 0; i < nbRowSad; ++i) {
        for (int j = 0; j < nbColSad; ++j) {
          dataType e = sadMatrix[i][j];
          if (e > thresholdSad && e < newThresholdSad) newThresholdSad = e;
        }
      }
      thresholdSad = newThresholdSad;

      // Apply Munkres again.
      for (int i = 0; i < nbRowSad; ++i) {
        for (int j = 0; j < nbColSad; ++j) {
          dataType e = sadMatrix[i][j];
          sadBottleneckMatrix[i][j] = (e > thresholdSad) ?
            std::numeric_limits<dataType>::max() : e;
        }
      }

      solverSad->setInput(nSad, nSad, (void*) sadBottleneckMatrix);
      solverSad->run<dataType>(sadMatchings);
      sadIter++;
    }

    if (sadBottleneckMatrix != NULL) {
      for (int r = 0; r < nSad; ++r)
        delete[] (sadBottleneckMatrix[r]);
      delete[] sadBottleneckMatrix;
      sadBottleneckMatrix = NULL;
    }

  }

  // Rebuild mappings.
  maxI = 0; minI = 0; sadI = 0;
  maxJ = 0; minJ = 0; sadJ = 0;

  vector<int> minMap1; vector<int> minMap2;
  vector<int> maxMap1; vector<int> maxMap2;
  vector<int> sadMap1; vector<int> sadMap2;

  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t1 = CTDiagram1->at(i);
    if (abs<dataType>(get<4>(t1)) < zeroThresh) continue;

    if (get<1>(t1) == NodeType::Local_minimum && get<3>(t1) == NodeType::Local_maximum) {
        maxMap1.push_back(i);
    } else {
      if (get<1>(t1) == NodeType::Local_minimum || get<3>(t1) == NodeType::Local_minimum)
        minMap1.push_back(i);
      if (get<1>(t1) == NodeType::Local_maximum || get<3>(t1) == NodeType::Local_maximum)
        maxMap1.push_back(i);
      if ((get<1>(t1) == NodeType::Saddle1 && get<3>(t1) == NodeType::Saddle2) ||
          (get<1>(t1) == NodeType::Saddle2 && get<3>(t1) == NodeType::Saddle1))
        sadMap1.push_back(i);
    }
  }

  for (int j = 0; j < d2Size; ++j) {
    diagramTuple t1 = CTDiagram2->at(j);
    if (abs<dataType>(get<4>(t1)) < zeroThresh) continue;

    if (get<1>(t1) == NodeType::Local_minimum && get<3>(t1) == NodeType::Local_maximum) {
      maxMap2.push_back(j);
    } else {
      if (get<1>(t1) == NodeType::Local_minimum || get<3>(t1) == NodeType::Local_minimum)
        minMap2.push_back(j);
      if (get<1>(t1) == NodeType::Local_maximum || get<3>(t1) == NodeType::Local_maximum)
        maxMap2.push_back(j);
      if ((get<1>(t1) == NodeType::Saddle1 && get<3>(t1) == NodeType::Saddle2) ||
          (get<1>(t1) == NodeType::Saddle2 && get<3>(t1) == NodeType::Saddle1))
        sadMap2.push_back(j);
    }
  }

  // Begin computation for unpaired vertices.
  dataType addedMinPersistence = 0;
  for (int i = 0, s = minMatchings->size(); i < s; ++i) {
    tuple<idVertex, idVertex, dataType> t = minMatchings->at(i);
    dataType val = abs<dataType>(get<2>(t));

    if (get<0>(t) >= (int) minMap1.size()) {
      addedMinPersistence = (wasserstein > 0 ?
        addedMinPersistence + std::pow(val, wasserstein) :
        std::max(val, addedMinPersistence));
    }
    else if (get<1>(t) >= (int) minMap2.size()) {
      addedMinPersistence = (wasserstein > 0 ?
        addedMinPersistence + std::pow(val, wasserstein) :
        std::max(val, addedMinPersistence));
    } else {
      matchings->push_back(make_tuple(minMap1[get<0>(t)], minMap2[get<1>(t)], val));
    }
  }

  dataType addedMaxPersistence = 0;
  for (int i = 0, s = maxMatchings->size(); i < s; ++i) {
    tuple<idVertex, idVertex, dataType> t = maxMatchings->at(i);
    dataType val = abs<dataType>(get<2>(t));

    if (get<0>(t) >= (int) maxMap1.size()) {
      addedMaxPersistence = (wasserstein > 0 ?
        addedMaxPersistence + std::pow(val, wasserstein) :
        std::max(val, addedMaxPersistence));
    }
    else if (get<1>(t) >= (int) maxMap2.size()) {
      addedMaxPersistence = (wasserstein > 0 ?
        addedMaxPersistence + std::pow(val, wasserstein) :
        std::max(val, addedMaxPersistence));
    } else {
      matchings->push_back(make_tuple(maxMap1[get<0>(t)], maxMap2[get<1>(t)], get<2>(t)));
    }
  }

  dataType addedSadPersistence = 0;
  for (int i = 0, s = sadMatchings->size(); i < s; ++i) {
    tuple<idVertex, idVertex, dataType> t = sadMatchings->at(i);
    dataType val = abs<dataType>(get<2>(t));

    if (get<0>(t) >= (int) sadMap1.size()) {
      addedSadPersistence = (wasserstein > 0 ?
        addedSadPersistence + std::pow(val, wasserstein) :
        std::max(val, addedSadPersistence));
    }
    else if (get<1>(t) >= (int) sadMap2.size()) {
      addedSadPersistence = (wasserstein > 0 ?
        addedSadPersistence + std::pow(val, wasserstein) :
        std::max(val, addedSadPersistence));
    } else {
      matchings->push_back(make_tuple(sadMap1[get<0>(t)], sadMap2[get<1>(t)], get<2>(t)));
    }
  }

  // Recompute matching weights for user-friendly distance.
  dataType d = 0;
  bool* paired1 = new bool[d1Size];
  bool* paired2 = new bool[d2Size];
  for (int b = 0; b < d1Size; ++b) paired1[b] = false;
  for (int b = 0; b < d2Size; ++b) paired2[b] = false;

  // cout << "[BottleneckDistance] Matchings (min/max): " << endl;
  // cout << endl << "[BottleneckDistance] Minima" << endl;
  // if (nbColMin < 59) solverMin->showCostMatrix<dataType>();
  // solverMin->showMaskMatrix();
  // cout << endl << "[BottleneckDistance] Maxima" << endl;
  // if (nbColMax < 50) solverMax->showCostMatrix<dataType>();
  // solverMax->showMaskMatrix();

  for (int m = 0, ms = (int) matchings->size(); m < ms; ++m)
  {
    tuple<idVertex, idVertex, dataType> t = matchings->at(m);
    int i = get<0>(t);
    int j = get<1>(t);

    diagramTuple t1 = CTDiagram1->at(i);
    diagramTuple t2 = CTDiagram2->at(j);
    dataType rX = get<6>(t1);
    dataType rY = get<10>(t1);
    dataType cX = get<6>(t2);
    dataType cY = get<10>(t2);
    dataType x = rX - cX;
    dataType y = rY - cY;

    paired1[i] = true;
    paired2[j] = true;
    dataType linfty = std::max(abs<dataType>(x), abs<dataType>(y));
    if (wasserstein > 0) {
      d += std::pow(linfty, wasserstein);
    } else {
      d = std::max(d, linfty);
    }
  }

  if (d1Size < d2Size) {
    for (int j = 0; j < d2Size; ++j) {
      if (paired2[j]) continue;
      diagramTuple t2 = CTDiagram2->at(j);
      dataType cDiff = abs<dataType>(get<4>(t2));
      d += wasserstein > 0 ? std::pow(cDiff, wasserstein) : 0;
    }
  } else if (d2Size < d1Size)  {
    for (int i = 0; i < d1Size; ++i) {
      if (paired1[i]) continue;
      diagramTuple t1 = CTDiagram1->at(i);
      dataType rDiff = abs<dataType>(get<4>(t1));
      d += wasserstein > 0 ? std::pow(rDiff, wasserstein) : 0;
      // Bottleneck => only maximum matching weight.
    }
  }

  d = wasserstein > 0 ?
      d + addedMaxPersistence + addedMinPersistence + addedSadPersistence :
      std::max(d, std::max(std::max(
          addedMaxPersistence, addedMinPersistence), addedSadPersistence));

  if (wasserstein > 0) {
    d = std::pow(d, (1.0 / (double) wasserstein));
  }

  stringstream msg;
  msg << "[BottleneckDistance] Computed distance " << d << endl;
  dMsg(std::cout, msg.str(), timeMsg);

  *distance = d;

  distance_ = (void*)(distance);
  return 0;
}

template <typename dataType>
bool BottleneckDistance::isValidMatching(
  const vector<tuple<idVertex, idVertex, dataType>>* matchings,
  dataType thresholdMin) const
{
  for (int i = 0, s = matchings->size(); i < s; ++i) {
    tuple<idVertex, idVertex, dataType> t = matchings->at(i);
    // Safe enough (matrix minima are substracted as a first
    // Munkres step, but std::numeric_limits should be high enough)
    if (get<2>(t) > thresholdMin)
      return false;
  }
  return true;
}

#endif // BOTTLENECKDISTANCE_H
