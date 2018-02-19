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
/// \sa ttkBottleneckDistance.cpp %for a usage example.

#ifndef _BOTTLENECKDISTANCE_H
#define _BOTTLENECKDISTANCE_H

#ifndef diagramTuple
#define diagramTuple tuple<ftm::idVertex, ftm::NodeType, ftm::idVertex, \
  ftm::NodeType, dataType, ftm::idVertex, \
  dataType, float, float, float, dataType, float, float, float>
#endif
#ifndef BNodeType
#define BNodeType ftm::NodeType
#define BLocalMax ftm::NodeType::Local_maximum
#define BLocalMin ftm::NodeType::Local_minimum
#define BSaddle1  ftm::NodeType::Saddle1
#define BSaddle2  ftm::NodeType::Saddle2
#define BIdVertex ftm::idVertex
#endif

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include                  <Munkres.h>

#include                  <string>
#include                  <tuple>
#include                  <functional>

namespace ttk {

  class BottleneckDistance : public Debug {

    public:
        
      BottleneckDistance():wasserstein_("inf") {};
      
      ~BottleneckDistance() {};

      template <typename dataType>
      int execute(bool usePersistenceMetric, double alpha);

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
      dataType getDistance() {
        return *static_cast<dataType*> (distance_);
      }

      template<typename type>
      static type abs(const type var) {
        return (var >= 0) ? var : -var;
      }

      template<typename type>
      static type abs_diff(const type var1, const type var2) {
        return (var1 > var2) ? var1 - var2 : var2 - var1;
      }
    
    protected:
    
      void                      *outputCT1_;
      void                      *outputCT2_;
      void                      *matchings_; // ids from CT1 to CT2
      void                      *distance_;

      string                    wasserstein_;
    
  private:

    template <typename dataType>
    int computeBottleneck(
      const vector<diagramTuple> *CTDiagram1,
      const vector<diagramTuple> *CTDiagram2,
      vector<matchingTuple> *matchings,
      bool usePersistenceMetric,
      double alpha);

    template <typename dataType>
    bool isValidMatching(
      const vector<matchingTuple>* matchings,
      dataType thresholdMin) const;

    template <typename dataType>
    double computeGeometricalRange(
      const vector<diagramTuple> *CTDiagram1,
      const vector<diagramTuple> *CTDiagram2,
      int d1Size,
      int d2Size) const;

    template <typename dataType>
    double computeMinimumRelevantPersistence(
      const vector<diagramTuple> *CTDiagram1,
      const vector<diagramTuple> *CTDiagram2,
      int d1Size,
      int d2Size) const;

    template <typename dataType>
    void computeMinMaxSaddleNumberAndMapping(
      const vector<diagramTuple> *CTDiagram,
      int dSize,
      int &nbMin,
      int &nbMax,
      int &nbSaddle,
      vector<int> *minMap,
      vector<int> *maxMap,
      vector<int> *sadMap,
      dataType zeroThresh);

    template <typename dataType>
    void buildCostMatrices(
      const vector<diagramTuple> *CTDiagram1,
      const vector<diagramTuple> *CTDiagram2,
      int d1Size,
      int d2Size,
      std::function<dataType (const diagramTuple, const diagramTuple)>& distanceFunction,
      std::function<dataType (const diagramTuple)>& diagonalDistanceFunction,
      double zeroThresh,
      dataType **minMatrix,
      dataType **maxMatrix,
      dataType **sadMatrix
    );

    template <typename dataType>
    dataType findInitialThreshold(
      int nbRow,
      int nbCol,
      const dataType **matrix
    ) const;

    template <typename dataType>
    void filterFromThreshold(
      dataType threshold,
      int nbRow,
      int nbCol,
      const dataType **matrix,
      dataType **bottleneckMatrix
    );

    template <typename dataType>
    void iterateSolving(
      vector<matchingTuple> *matchings,
      dataType threshold,
      int nbRow,
      int nbCol,
      const dataType **matrix,
      dataType **bottleneckMatrix,
      Munkres *solver);

    template <typename dataType>
    void solvePWasserstein(
      int nbRow,
      int nbCol,
      dataType **matrix,
      vector<matchingTuple> *matchings,
      Munkres *solver);

    template <typename dataType>
    void solveInfinityWasserstein(
      int nbRow,
      int nbCol,
      dataType **matrix,
      vector<matchingTuple> *matchings,
      Munkres *solver);

    template <typename dataType>
    dataType buildMappings(
      vector<matchingTuple> inputMatchings,
      vector<matchingTuple> *outputMatchings,
      vector<int> map1,
      vector<int> map2,
      int wasserstein);

  };

}

#include <BottleneckDistanceImpl.h>
#include <BottleneckDistanceMainImpl.h>

#endif // _H
