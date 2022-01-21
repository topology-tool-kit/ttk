/// \ingroup base
/// \class ttk::ReebSpace
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2015.
///
/// \brief TTK processing package that efficiently computes the Reeb space of
/// bivariate volumetric data.
///
/// The Reeb space is a useful topological abstraction of bivariate scalar
/// fields for data segmentation purposes. Intuitively, it allows the automatic
/// separation of volumetric regions that project to the same areas in the
/// range. This class also implements a simplification heuristic for progressive
/// coarsening. Used in conjunction with continuous scatterplots, this class
/// enables continuous scattterplot peeling for instance.
///
/// \param dataTypeU Data type of the input first component field (char, float,
/// etc.)
/// \param dataTypeV Data type of the input second component field (char, float,
/// etc.)
///
/// \b Related \b publication \n
/// "Jacobi Fiber Surfaces for Bivariate Reeb Space Computation" \n
/// Julien Tierny, Hamish Carr \n
/// Proc. of IEEE VIS 2016.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa ttkReebSpace.cpp %for a usage example.

#pragma once

// base code includes
#include <FiberSurface.h>
#include <Geometry.h>
#include <JacobiSet.h>
#include <Triangulation.h>

#include <map>
#include <set>

namespace ttk {

  class ReebSpace : virtual public Debug {
  public:
    enum class SimplificationCriterion {
      domainVolume, // 0
      rangeArea, // 1
      hyperVolume // 2
    };

    struct Sheet0 {
      SimplexId vertexId_{}, pruned_{};
      // 0: extrema-sheet, 1: saddle-sheet
      char type_{};
      std::vector<SimplexId> sheet1List_{};
      std::vector<SimplexId> sheet3List_{};
    };

    struct Sheet1 {
      bool hasSaddleEdges_{}, pruned_{};
      std::vector<SimplexId> edgeList_{};
      // 0: extrema-sheet, 1: saddle-sheet (can be inferred by the number
      // of 2-sheets attached to it?)
      std::vector<SimplexId> sheet0List_{};
      // NB: the corresponding 2-sheet should have the same global
      // indentifier.
      std::vector<SimplexId> sheet3List_{};
    };

    struct Sheet2 {
      bool pruned_{};
      SimplexId sheet1Id_{};
      // for each point, x, y, z coordinates
      // this is how coincident point will be merged afterwards
      // this list is meant to be temporary. after the merge, we'll use a
      // global list
      std::vector<std::vector<FiberSurface::Vertex>> vertexList_{};
      std::vector<std::vector<FiberSurface::Triangle>> triangleList_{};
      std::vector<SimplexId> sheet3List_{};
    };

    struct Sheet3 {
      SimplexId Id_{}, simplificationId_{}, preMerger_{};
      bool pruned_{};
      double domainVolume_{}, rangeArea_{}, hyperVolume_{};
      std::vector<SimplexId> vertexList_{};
      std::vector<SimplexId> tetList_{};
      std::vector<SimplexId> sheet0List_{};
      std::vector<SimplexId> sheet1List_{};
      std::vector<SimplexId> sheet2List_{};
      std::vector<SimplexId> sheet3List_{};
      std::vector<SimplexId> preMergedSheets_{};
    };

    ReebSpace();

    inline bool empty() const {
      return (currentData_.vertex2sheet0_.size() == 0);
    }

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int execute(const dataTypeU *const uField,
                       const dataTypeV *const vField,
                       const triangulationType &triangulation);

    inline const Sheet0 *get0sheet(const SimplexId &sheetId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((sheetId < 0)
         || (sheetId >= (SimplexId)currentData_.sheet0List_.size()))
        return NULL;
#endif

      return &(currentData_.sheet0List_[sheetId]);
    }

    inline const Sheet1 *get1sheet(const SimplexId &sheetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((sheetId < 0)
         || (sheetId >= (SimplexId)currentData_.sheet1List_.size()))
        return NULL;
#endif

      return &(currentData_.sheet1List_[sheetId]);
    }

    // warning, these are the originals
    inline const Sheet2 *get2sheet(const SimplexId &sheetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((sheetId < 0)
         || (sheetId >= (SimplexId)originalData_.sheet2List_.size()))
        return NULL;
#endif

      return &(originalData_.sheet2List_[sheetId]);
    }

    inline const Sheet3 *get3sheet(const SimplexId &sheetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((sheetId < 0)
         || (sheetId >= (SimplexId)currentData_.sheet3List_.size()))
        return NULL;
#endif

      return &(currentData_.sheet3List_[sheetId]);
    }

    inline const std::vector<SimplexId> *get0sheetSegmentation() const {
      return &currentData_.vertex2sheet0_;
    }

    const std::vector<SimplexId> *get1sheetSegmentation() const {
      return &currentData_.edge2sheet1_;
    }

    const std::vector<SimplexId> *get3sheetVertexSegmentation() const {
      return &currentData_.vertex2sheet3_;
    }

    const std::vector<SimplexId> *get3sheetTetSegmentation() const {
      return &currentData_.tet2sheet3_;
    }

    const std::vector<char> *getEdgeTypes() const {
      return &currentData_.edgeTypes_;
    }

    inline const std::vector<FiberSurface::Vertex> *
      getFiberSurfaceVertices() const {
      return &(fiberSurfaceVertexList_);
    }

    inline int getJacobi2Edge(const SimplexId &jacobiEdgeId) const {
      if((jacobiEdgeId < 0)
         || (jacobiEdgeId >= (SimplexId)jacobi2edges_.size()))
        return -1;
      return jacobi2edges_[jacobiEdgeId];
    }

    inline SimplexId getNumberOf2sheets() const {
      return currentData_.sheet2List_.size();
    }

    //       inline std::vector<long long int>* getSheetTriangulationCells(){
    //         return &sheet3cells_;
    //       }
    //
    //       inline std::vector<float>* getSheetTriangulationPoints(){
    //         return &sheet3points_;
    //       }

    template <class dataTypeU, class dataTypeV>
    inline int
      perturbate(const dataTypeU *const uField,
                 const dataTypeV *const vField,
                 const dataTypeU &uEpsilon = Geometry::powIntTen(-DBL_DIG),
                 const dataTypeV &vEpsilon = Geometry::powIntTen(-DBL_DIG));

    inline void setExpand3Sheets(const bool &onOff) {
      expand3sheets_ = onOff;
    }

    inline bool setRangeDrivenOctree(const bool &onOff) {
      if(onOff != withRangeDrivenOctree_) {
        withRangeDrivenOctree_ = onOff;
        return true;
      }

      withRangeDrivenOctree_ = onOff;
      return false;
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p sosOffsetsU buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setSosOffsetsU(const SimplexId *const sosOffsetsU) {
      sosOffsetsU_ = sosOffsetsU;
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p sosOffsetsV buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setSosOffsetsV(const SimplexId *const sosOffsetsV) {
      sosOffsetsV_ = sosOffsetsV;
    }

    inline void setTetNumber(const SimplexId &tetNumber) {
      tetNumber_ = tetNumber;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    // WARNING: if you plan to use the range driven octree, make sure
    // that you provided pointers to the u and v fields.
    inline int
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation) {
        triangulation->preconditionVertexStars();
        triangulation->preconditionEdges();
        triangulation->preconditionVertexEdges();

        // trigger the jacobiSet pre-processing on the triangulation.
        jacobiSet_.preconditionTriangulation(triangulation);

        // trigger the fiberSurface pre-processing on the triangulation.
        fiberSurface_.preconditionTriangulation(triangulation);

        vertexNumber_ = triangulation->getNumberOfVertices();
        edgeNumber_ = triangulation->getNumberOfEdges();
        tetNumber_ = triangulation->getNumberOfCells();
      }

      return 0;
    }

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int simplify(const dataTypeU *const uField,
                        const dataTypeV *const vField,
                        const triangulationType &triangulation,
                        const double &simplificationThreshold,
                        const SimplificationCriterion &criterion
                        = SimplificationCriterion::rangeArea);

  private:
    // output segmentation
    struct ReebSpaceData {
      SimplificationCriterion simplificationCriterion_{};
      double simplificationThreshold_{};

      std::vector<SimplexId> edge2sheet1_{};
      std::vector<char> edgeTypes_{};
      std::vector<SimplexId> tet2sheet3_{};
      std::vector<SimplexId> vertex2sheet0_{};
      std::vector<SimplexId> vertex2sheet3_{};

      // structure
      std::vector<Sheet0> sheet0List_{};
      std::vector<Sheet1> sheet1List_{};
      std::vector<Sheet2> sheet2List_{};
      std::vector<Sheet3> sheet3List_{};

      //         std::vector<float>         sheet3points_;
      //         std::vector<long long int> sheet3cells_;
    };

  protected:
    template <typename triangulationType>
    int compute1sheetsOnly(
      const std::vector<std::pair<SimplexId, char>> &jacobiSet,
      std::vector<std::pair<SimplexId, SimplexId>> &jacobiSetClassification,
      const triangulationType &triangulation);

    template <typename triangulationType>
    int compute1sheets(
      const std::vector<std::pair<SimplexId, char>> &jacobiSet,
      std::vector<std::pair<SimplexId, SimplexId>> &jacobiSetClassification,
      const triangulationType &triangulation);

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int compute2sheets(
      const std::vector<std::pair<SimplexId, SimplexId>> &jacobiEdges,
      const dataTypeU *const uField,
      const dataTypeV *const vField,
      const triangulationType &triangulation);

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int compute2sheetChambers(const dataTypeU *const uField,
                                     const dataTypeV *const vField,
                                     const triangulationType &triangulation);

    template <typename triangulationType>
    int compute3sheet(
      const SimplexId &vertexId,
      const std::vector<std::vector<std::vector<SimplexId>>> &tetTriangles,
      const triangulationType &triangulation);

    template <typename triangulationType>
    int compute3sheets(
      std::vector<std::vector<std::vector<SimplexId>>> &tetTriangles,
      const triangulationType &triangulation);

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int
      computeGeometricalMeasures(Sheet3 &sheet,
                                 const dataTypeU *const uField,
                                 const dataTypeV *const vField,
                                 const triangulationType &triangulation);

    int connect3sheetTo0sheet(ReebSpaceData &data,
                              const SimplexId &sheet3Id,
                              const SimplexId &sheet0Id);

    int connect3sheetTo1sheet(ReebSpaceData &data,
                              const SimplexId &sheet3Id,
                              const SimplexId &sheet1Id);

    int connect3sheetTo2sheet(ReebSpaceData &data,
                              const SimplexId &sheet3Id,
                              const SimplexId &sheet2Id);

    int connect3sheetTo3sheet(ReebSpaceData &data,
                              const SimplexId &sheet3Id,
                              const SimplexId &otherSheet3Id);

    template <typename triangulationType>
    int connectSheets(const triangulationType &triangulation);

    int disconnect1sheetFrom0sheet(ReebSpaceData &data,
                                   const SimplexId &sheet1Id,
                                   const SimplexId &sheet0Id,
                                   const SimplexId &biggerId);

    int disconnect3sheetFrom0sheet(ReebSpaceData &data,
                                   const SimplexId &sheet3Id,
                                   const SimplexId &sheet0Id);

    template <typename triangulationType>
    int disconnect3sheetFrom1sheet(ReebSpaceData &data,
                                   const SimplexId &sheet3Id,
                                   const SimplexId &sheet1Id,
                                   const SimplexId &biggerId,
                                   const triangulationType &triangulation);

    int disconnect3sheetFrom2sheet(ReebSpaceData &data,
                                   const SimplexId &sheet3Id,
                                   const SimplexId &sheet2Id);

    int disconnect3sheetFrom3sheet(ReebSpaceData &data,
                                   const SimplexId &sheet3Id,
                                   const SimplexId &other3SheetId);

    template <typename triangulationType>
    int flush(const triangulationType &triangulation);

    template <typename triangulationType>
    int mergeSheets(const SimplexId &smallerId,
                    const SimplexId &biggerId,
                    const triangulationType &triangulation);

    int preMergeSheets(const SimplexId &sheetId0, const SimplexId &sheetId1);

    int prepareSimplification();

    int printConnectivity(const ReebSpaceData &data) const;

    template <typename triangulationType>
    int simplifySheets(const double &simplificationThreshold,
                       const SimplificationCriterion &simplificationCriterion,
                       const triangulationType &triangulation);

    template <typename triangulationType>
    int simplifySheet(const SimplexId &sheetId,
                      const SimplificationCriterion &simplificationCriterion,
                      const triangulationType &triangulation);

    //       int triangulateTetrahedron(const int &tetId,
    //         const std::vector<std::vector<int> > &triangles,
    //         std::vector<long long int> &outputTets);
    //
    //       int triangulateThreeSheets();

    SimplexId vertexNumber_{0}, edgeNumber_{0}, tetNumber_{0};
    double totalArea_{-1}, totalVolume_{-1}, totalHyperVolume_{-1};

    const SimplexId *sosOffsetsU_{}, *sosOffsetsV_{};

    bool hasConnectedSheets_{false}, expand3sheets_{true},
      withRangeDrivenOctree_{true};
    ReebSpaceData originalData_{}, currentData_{};

    // information that does not get simplified
    std::vector<std::pair<SimplexId, char>> jacobiSetEdges_{};
    std::vector<SimplexId> jacobi2edges_{};

    FiberSurface fiberSurface_{};
    JacobiSet jacobiSet_{};
    std::vector<FiberSurface::Vertex> fiberSurfaceVertexList_{};
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <ReebSpace.cpp>

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::ReebSpace::execute(const dataTypeU *const uField,
                                   const dataTypeV *const vField,
                                   const triangulationType &triangulation) {

  flush(triangulation);
  fiberSurface_.setInputField(uField, vField);

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  fiberSurface_.flushOctree();
  if(withRangeDrivenOctree_) {
    fiberSurface_.buildOctree<dataTypeU, dataTypeV>(&triangulation);
  }
#endif

  Timer t;

  // 1) compute the jacobi set
  jacobiSet_.setSosOffsetsU(sosOffsetsU_);
  jacobiSet_.setSosOffsetsV(sosOffsetsV_);
  jacobiSet_.execute(jacobiSetEdges_, uField, vField, triangulation);

  // 2) compute the list saddle 1-sheets
  // + list of saddle 0-sheets
  std::vector<std::pair<SimplexId, SimplexId>> jacobiSetClassification;
  compute1sheetsOnly(jacobiSetEdges_, jacobiSetClassification, triangulation);
  // at this stage, jacobiSetClassification contains the list of saddle edges
  // along with their 1-sheet Id.

  compute2sheets(jacobiSetClassification, uField, vField, triangulation);
  //   compute2sheetChambers<dataTypeU, dataTypeV>();

  std::vector<std::vector<std::vector<SimplexId>>> tetTriangles;
  compute3sheets(tetTriangles, triangulation);

  this->printMsg(
    "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);

  // post-process for further interaction
  if((totalArea_ == -1) || (totalVolume_ == -1) || (totalHyperVolume_ == -1)) {

    Timer tm;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
      computeGeometricalMeasures(
        originalData_.sheet3List_[i], uField, vField, triangulation);
    }

    for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
      totalArea_ += originalData_.sheet3List_[i].rangeArea_;
      totalVolume_ += originalData_.sheet3List_[i].domainVolume_;
      totalHyperVolume_ += originalData_.sheet3List_[i].hyperVolume_;
    }

    this->printMsg("Computed geometrical measures", 1.0, tm.getElapsedTime(),
                   this->threadNumber_);
  }

  fiberSurface_.finalize<dataTypeU, dataTypeV>();

  prepareSimplification();

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::ReebSpace::compute2sheets(
  const std::vector<std::pair<SimplexId, SimplexId>> &jacobiEdges,
  const dataTypeU *const uField,
  const dataTypeV *const vField,
  const triangulationType &triangulation) {

  Timer t;

  // at this point, they have exactly the same size
  originalData_.sheet2List_.resize(originalData_.sheet1List_.size());
  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {
    originalData_.sheet2List_[i].sheet1Id_ = i;
    originalData_.sheet2List_[i].pruned_ = false;
    originalData_.sheet2List_[i].vertexList_.resize(
      originalData_.sheet1List_[originalData_.sheet2List_[i].sheet1Id_]
        .edgeList_.size());
    originalData_.sheet2List_[i].triangleList_.resize(
      originalData_.sheet1List_[originalData_.sheet2List_[i].sheet1Id_]
        .edgeList_.size());
    for(size_t j = 0; j < originalData_.sheet2List_[i].vertexList_.size();
        j++) {
      originalData_.sheet2List_[i].vertexList_[j].clear();
      originalData_.sheet2List_[i].triangleList_[j].clear();
    }
  }

  fiberSurface_.setGlobalVertexList(&fiberSurfaceVertexList_);
  fiberSurface_.setPolygonEdgeNumber(jacobiEdges.size());

  std::vector<SimplexId> edge2polygonEdgeId(edgeNumber_, -1);
  jacobi2edges_.resize(jacobiEdges.size());

  for(size_t i = 0; i < jacobiEdges.size(); i++) {

    SimplexId edgeId = jacobiEdges[i].first;

    edge2polygonEdgeId[edgeId] = i;
    jacobi2edges_[i] = edgeId;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {

    for(size_t j = 0;
        j < originalData_.sheet1List_[originalData_.sheet2List_[i].sheet1Id_]
              .edgeList_.size();
        j++) {

      SimplexId edgeId
        = originalData_.sheet1List_[originalData_.sheet2List_[i].sheet1Id_]
            .edgeList_[j];

      fiberSurface_.setTriangleList(
        edge2polygonEdgeId[edgeId],
        &(originalData_.sheet2List_[i].triangleList_[j]));
      fiberSurface_.setVertexList(
        edge2polygonEdgeId[edgeId],
        &(originalData_.sheet2List_[i].vertexList_[j]));
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < jacobiEdges.size(); i++) {

    SimplexId edgeId = jacobiEdges[i].first;

    std::pair<double, double> rangePoint0, rangePoint1;

    SimplexId vertexId0 = -1, vertexId1 = -1;
    triangulation.getEdgeVertex(edgeId, 0, vertexId0);
    triangulation.getEdgeVertex(edgeId, 1, vertexId1);

    rangePoint0.first = uField[vertexId0];
    rangePoint0.second = vField[vertexId0];

    rangePoint1.first = uField[vertexId1];
    rangePoint1.second = vField[vertexId1];

    if(originalData_.edgeTypes_[edgeId] == 1) {
      std::vector<SimplexId> edgeSeeds(
        triangulation.getEdgeStarNumber(edgeId), -1);
      for(size_t j = 0; j < edgeSeeds.size(); j++) {
        triangulation.getEdgeStar(edgeId, j, edgeSeeds[j]);
      }

      fiberSurface_.computeContour<dataTypeU, dataTypeV>(
        rangePoint0, rangePoint1, edgeSeeds, &triangulation,
        edge2polygonEdgeId[edgeId]);
    } else {
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
      if(withRangeDrivenOctree_) {
        fiberSurface_.computeSurfaceWithOctree<dataTypeU, dataTypeV>(
          rangePoint0, rangePoint1, &triangulation, edge2polygonEdgeId[edgeId]);
      } else {
        fiberSurface_.computeSurface<dataTypeU, dataTypeV>(
          rangePoint0, rangePoint1, &triangulation, edge2polygonEdgeId[edgeId]);
      }
#else
      fiberSurface_.computeSurface<dataTypeU, dataTypeV>(
        rangePoint0, rangePoint1, edge2polygonEdgeId[edgeId]);
#endif
    }
  }

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  //   for(int i = 0; i < (int) originalData_.sheet2List_.size(); i++){
  //
  //     int sheet1Id = originalData_.sheet2List_[i].sheet1Id_;
  //
  //     std::vector<bool> inList(tetNumber_, false);
  //     std::vector<int> seedList;
  //     std::vector<std::pair<std::pair<double, double>, std::pair<double,
  //      double> > > edgeList;
  //     std::vector<int> jacobiEdgeIdList;
  //
  //     for(int j = 0;
  //       j < (int) originalData_.sheet1List_[sheet1Id].edgeList_.size(); j++){
  //       int edgeId = originalData_.sheet1List_[sheet1Id].edgeList_[j];
  //
  //       if(originalData_.edgeTypes_[edgeId] == 1){
  //         for(int k = 0; k < (int) (*edgeStars_)[edgeId].size(); k++){
  //           int tetId = (*edgeStars_)[edgeId][k];
  //           if(!inList[tetId]){
  //             seedList.push_back(tetId);
  //           }
  //         }
  //       }
  //
  //       if((originalData_.edgeTypes_[edgeId] != 1)
  //         &&(originalData_.edgeTypes_[edgeId] != -1)){
  //         edgeList.resize(edgeList.size() + 1);
  //
  //         edgeList.back().first.first =
  //           ((dataTypeU *) uField_)[(*edgeList_)[edgeId].first];
  //         edgeList.back().first.second =
  //           ((dataTypeV *) vField_)[(*edgeList_)[edgeId].first];
  //
  //         edgeList.back().second.first =
  //           ((dataTypeU *) uField_)[(*edgeList_)[edgeId].second];
  //         edgeList.back().second.second =
  //           ((dataTypeV *) vField_)[(*edgeList_)[edgeId].second];
  //
  //         jacobiEdgeIdList.push_back(edge2polygonEdgeId[edgeId]);
  //       }
  //     }
  //
  //     if(seedList.size()){
  //       fiberSurface_.computeContour<dataTypeU, dataTypeV>(
  //         edgeList, seedList,
  //         &jacobiEdgeIdList);
  //     }
  //   }

  this->printMsg(
    "Computed fiber surfaces", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::ReebSpace::compute2sheetChambers(
  const dataTypeU *const uField,
  const dataTypeV *const vField,
  const triangulationType &triangulation) {

  this->printWrn("Computing chambers' pre-images.");
  this->printWrn("This will take a LONG time.");

  Timer t;

  // at this point, they have exactly the same size
  originalData_.sheet2List_.resize(threadNumber_);
  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {
    originalData_.sheet2List_[i].sheet1Id_ = i;
    originalData_.sheet2List_[i].pruned_ = false;
    originalData_.sheet2List_[i].vertexList_.resize(1);
    originalData_.sheet2List_[i].triangleList_.resize(1);
  }

  fiberSurface_.setGlobalVertexList(&fiberSurfaceVertexList_);
  fiberSurface_.setTetNumber(tetNumber_);
  fiberSurface_.setPolygonEdgeNumber(threadNumber_);

  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {

    fiberSurface_.setTriangleList(
      i, &(originalData_.sheet2List_[i].triangleList_[0]));
    fiberSurface_.setVertexList(
      i, &(originalData_.sheet2List_[i].vertexList_[0]));
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)edgeNumber_; i++) {

#ifdef TTK_ENABLE_OPENMP
    ThreadId threadId = omp_get_thread_num();
#else
    ThreadId threadId = 0;
#endif

    SimplexId edgeId = i;

    std::pair<double, double> rangePoint0, rangePoint1;

    SimplexId vertexId0 = -1, vertexId1 = -1;
    triangulation.getEdgeVertex(edgeId, 0, vertexId0);
    triangulation.getEdgeVertex(edgeId, 1, vertexId1);

    rangePoint0.first = uField[vertexId0];
    rangePoint0.second = vField[vertexId0];

    rangePoint1.first = uField[vertexId1];
    rangePoint1.second = vField[vertexId1];

    fiberSurface_.computeSurface<dataTypeU, dataTypeV>(
      rangePoint0, rangePoint1, threadId);

    // clear the memory now.. otherwise we'll swap and never end
    originalData_.sheet2List_[threadId].triangleList_[0].clear();
    originalData_.sheet2List_[threadId].vertexList_[0].clear();
  }

  this->printMsg("Computed chambers pre-image boundaries", 1.0,
                 t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::ReebSpace::computeGeometricalMeasures(
  Sheet3 &sheet,
  const dataTypeU *const uField,
  const dataTypeV *const vField,
  const triangulationType &triangulation) {

  sheet.domainVolume_ = 0;
  sheet.rangeArea_ = 0;
  sheet.hyperVolume_ = 0;

  for(size_t i = 0; i < sheet.tetList_.size(); i++) {

    SimplexId tetId = sheet.tetList_[i];
    std::vector<std::pair<double, double>> domainBox, rangeBox;
    std::vector<std::vector<float>> domainPoints(4), rangePoints(4);

    for(int j = 0; j < 4; j++) {
      domainPoints[j].resize(3);
      rangePoints[j].resize(2);

      SimplexId vertexId = -1;
      triangulation.getCellVertex(tetId, j, vertexId);

      triangulation.getVertexPoint(
        vertexId, domainPoints[j][0], domainPoints[j][1], domainPoints[j][2]);

      rangePoints[j][0] = uField[vertexId];
      rangePoints[j][1] = vField[vertexId];
    }

    Geometry::getBoundingBox(domainPoints, domainBox);
    Geometry::getBoundingBox(rangePoints, rangeBox);

    sheet.domainVolume_ += (domainBox[0].second - domainBox[0].first)
                           * (domainBox[1].second - domainBox[1].first)
                           * (domainBox[2].second - domainBox[2].first);

    sheet.rangeArea_ += (rangeBox[0].second - rangeBox[0].first)
                        * (rangeBox[1].second - rangeBox[1].first);
  }

  if(sheet.domainVolume_) {
    sheet.hyperVolume_ = sheet.rangeArea_ / sheet.domainVolume_;
  } else {
    sheet.hyperVolume_ = 0;
  }

  return 0;
}

template <class dataTypeU, class dataTypeV>
inline int ttk::ReebSpace::perturbate(const dataTypeU *const uField,
                                      const dataTypeV *const vField,
                                      const dataTypeU &uEpsilon,
                                      const dataTypeV &vEpsilon) {

  jacobiSet_.setVertexNumber(vertexNumber_);
  jacobiSet_.perturbate(uField, vField, uEpsilon, vEpsilon);

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::ReebSpace::simplify(const dataTypeU *const uField,
                                    const dataTypeV *const vField,
                                    const triangulationType &triangulation,
                                    const double &simplificationThreshold,
                                    const SimplificationCriterion &criterion) {

  if((totalArea_ == -1) || (totalVolume_ == -1) || (totalHyperVolume_ == -1)) {

    Timer t;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
      computeGeometricalMeasures(
        originalData_.sheet3List_[i], uField, vField, triangulation);
    }

    for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
      totalArea_ += originalData_.sheet3List_[i].rangeArea_;
      totalVolume_ += originalData_.sheet3List_[i].domainVolume_;
      totalHyperVolume_ += originalData_.sheet3List_[i].hyperVolume_;
    }

    this->printMsg("Computed geometrical measures", 1.0, t.getElapsedTime(),
                   this->threadNumber_);
  }

  if(!hasConnectedSheets_) {
    connectSheets(triangulation);
    prepareSimplification();
  }

  {
    std::stringstream msg;
    msg << "Simplifying (";
    switch(criterion) {
      case SimplificationCriterion::domainVolume:
        msg << "'Domain Volume'";
        break;
      case SimplificationCriterion::rangeArea:
        msg << "'Range Area'";
        break;
      case SimplificationCriterion::hyperVolume:
        msg << "'HyperVolume'";
        break;
    }
    msg << ", thr: " << simplificationThreshold << ").";
    this->printMsg(msg.str());
  }

  if(!((criterion == currentData_.simplificationCriterion_)
       && (simplificationThreshold > currentData_.simplificationThreshold_))) {
    prepareSimplification();
  }

  simplifySheets(simplificationThreshold, criterion, triangulation);

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::compute1sheetsOnly(
  const std::vector<std::pair<SimplexId, char>> &jacobiSet,
  std::vector<std::pair<SimplexId, SimplexId>> &jacobiClassification,
  const triangulationType &triangulation) {

  // bfs on the saddle jacobi edges only to identify saddle 1-sheets as well as
  // saddle 0-sheets

  Timer t;

  // alloc
  jacobiClassification.reserve(jacobiSet.size());

  // markup the saddle jacobi and visisted edges
  for(size_t i = 0; i < jacobiSet.size(); i++) {
    originalData_.edgeTypes_[jacobiSet[i].first] = jacobiSet[i].second;
  }

  std::vector<bool> visitedEdges(triangulation.getNumberOfEdges(), false);

  for(size_t i = 0; i < jacobiSet.size(); i++) {

    if(visitedEdges[jacobiSet[i].first] == false) {

      SimplexId sheet1Id = originalData_.sheet1List_.size();
      originalData_.sheet1List_.resize(originalData_.sheet1List_.size() + 1);
      originalData_.sheet1List_.back().hasSaddleEdges_ = false;
      originalData_.sheet1List_.back().pruned_ = false;

      std::queue<SimplexId> edgeQueue;
      edgeQueue.push(jacobiSet[i].first);

      do {

        SimplexId edgeId = edgeQueue.front();
        edgeQueue.pop();

        if(!visitedEdges[edgeId]) {

          jacobiClassification.push_back(
            std::pair<SimplexId, SimplexId>(edgeId, sheet1Id));

          originalData_.sheet1List_.back().edgeList_.push_back(edgeId);
          originalData_.edge2sheet1_[edgeId] = sheet1Id;

          if(originalData_.edgeTypes_[edgeId] == 1) {
            originalData_.sheet1List_.back().hasSaddleEdges_ = true;
          }

          SimplexId vertexId0 = -1;
          triangulation.getEdgeVertex(edgeId, 0, vertexId0);
          SimplexId vertexId1 = -1;
          triangulation.getEdgeVertex(edgeId, 1, vertexId1);

          SimplexId neighborNumber = 0;

          SimplexId vertexEdgeNumber
            = triangulation.getVertexEdgeNumber(vertexId0);

          for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

            SimplexId vertexEdgeId = -1;
            triangulation.getVertexEdge(vertexId0, j, vertexEdgeId);

            if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
              if(vertexEdgeId != edgeId) {

                neighborNumber++;

                if(!visitedEdges[vertexEdgeId]) {
                  edgeQueue.push(vertexEdgeId);
                }
              }
            }
          }
          if((neighborNumber > 2)
             && (originalData_.vertex2sheet0_[vertexId0] == -1)) {
            // vertexId0 is a 0-sheet
            // mark up the segmentation
            originalData_.vertex2sheet0_[vertexId0]
              = originalData_.sheet0List_.size();

            originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                             + 1);
            originalData_.sheet0List_.back().vertexId_ = vertexId0;
            originalData_.sheet0List_.back().type_ = 1;
            originalData_.sheet0List_.back().pruned_ = false;
          }

          neighborNumber = 0;
          vertexEdgeNumber = triangulation.getVertexEdgeNumber(vertexId1);
          for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

            SimplexId vertexEdgeId = -1;
            triangulation.getVertexEdge(vertexId1, j, vertexEdgeId);

            if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
              if(vertexEdgeId != edgeId) {

                neighborNumber++;

                if(!visitedEdges[vertexEdgeId]) {
                  edgeQueue.push(vertexEdgeId);
                }
              }
            }
          }

          if((neighborNumber > 2)
             && (originalData_.vertex2sheet0_[vertexId1] == -1)) {
            // vertexId1 is a 0-sheet
            // mark up the segmentation

            originalData_.vertex2sheet0_[vertexId1]
              = originalData_.sheet0List_.size();

            originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                             + 1);
            originalData_.sheet0List_.back().vertexId_ = vertexId1;
            originalData_.sheet0List_.back().type_ = 1;
            originalData_.sheet0List_.back().pruned_ = false;
          }
        }

        visitedEdges[edgeId] = true;

      } while(edgeQueue.size());
    }
  }

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#1-sheets", std::to_string(originalData_.sheet1List_.size())},
    {"#0-sheets", std::to_string(originalData_.sheet0List_.size())}});
  this->printMsg(
    "Extracted 0- and 1-sheets", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::compute1sheets(
  const std::vector<std::pair<SimplexId, char>> &jacobiSet,
  std::vector<std::pair<SimplexId, SimplexId>> &jacobiClassification,
  const triangulationType &triangulation) {

  // bfs on the saddle jacobi edges only to identify saddle 1-sheets as well as
  // saddle 0-sheets

  Timer t;

  // alloc
  jacobiClassification.reserve(jacobiSet.size());

  // markup the saddle jacobi and visisted edges
  for(size_t i = 0; i < jacobiSet.size(); i++) {
    originalData_.edgeTypes_[jacobiSet[i].first] = jacobiSet[i].second;
  }

  std::vector<bool> visitedEdges(triangulation.getNumberOfEdges(), false);
  for(size_t i = 0; i < jacobiSet.size(); i++) {

    if(/*(saddleEdge[jacobiSet[i].first] == 1)&&*/
       (visitedEdges[jacobiSet[i].first] == false)) {
      // saddle, non-visited edge

      // we have a seed here
      SimplexId sheet1Id = originalData_.sheet1List_.size();
      originalData_.sheet1List_.resize(originalData_.sheet1List_.size() + 1);
      originalData_.sheet1List_.back().hasSaddleEdges_ = false;
      originalData_.sheet1List_.back().pruned_ = false;

      std::queue<SimplexId> edgeQueue;
      edgeQueue.push(jacobiSet[i].first);

      do {

        SimplexId edgeId = edgeQueue.front();
        edgeQueue.pop();

        if(!visitedEdges[edgeId]) {

          jacobiClassification.push_back(
            std::pair<SimplexId, SimplexId>(edgeId, sheet1Id));

          if(originalData_.edgeTypes_[edgeId] == 1) {
            // saddle edge
            originalData_.sheet1List_.back().hasSaddleEdges_ = true;
          }
          originalData_.sheet1List_.back().edgeList_.push_back(edgeId);
          originalData_.edge2sheet1_[edgeId] = sheet1Id;

          // next grab its neighbors
          // if saddleEdge and not visited continue (only if one)
          SimplexId vertexId0 = -1;
          triangulation.getEdgeVertex(edgeId, 0, vertexId0);
          SimplexId vertexId1 = -1;
          triangulation.getEdgeVertex(edgeId, 1, vertexId1);

          if(originalData_.vertex2sheet0_[vertexId0]
             >= static_cast<SimplexId>(originalData_.sheet0List_.size())) {
            // WEIRD BUG after multiple runs
            originalData_.vertex2sheet0_[vertexId0] = -1;
          }

          // make sure these are not already visited 0-sheets
          if(originalData_.vertex2sheet0_[vertexId0] == -1) {

            // inspect neighbor edges
            SimplexId neighborNumber = 0;
            SimplexId nextEdgeId = -1;
            SimplexId vertexEdgeNumber
              = triangulation.getVertexEdgeNumber(vertexId0);
            for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

              SimplexId vertexEdgeId = -1;
              triangulation.getVertexEdge(vertexId0, j, vertexEdgeId);

              if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
                neighborNumber++;
                if(vertexEdgeId != edgeId
                   /*&&(saddleEdge[(*vertexEdgeList_)[vertexId0][j]] == 1)*/) {
                  nextEdgeId = vertexEdgeId;
                }
              }
            }

            if((neighborNumber == 2) && (nextEdgeId != -1)) {
              // this is a regular vertex and we can add the next edge to the
              // queue if it's not been visited before.
              if(!visitedEdges[nextEdgeId])
                edgeQueue.push(nextEdgeId);
            } else {
              // we hit a 0-sheet that hasn't been visited before.
              // let's create it

              // mark up the segmentation
              originalData_.vertex2sheet0_[vertexId0]
                = originalData_.sheet0List_.size();

              originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                               + 1);
              originalData_.sheet0List_.back().vertexId_ = vertexId0;
              originalData_.sheet0List_.back().type_ = 1;
              originalData_.sheet0List_.back().pruned_ = false;
            }
          }

          if(originalData_.vertex2sheet0_[vertexId0] != -1) {
            // we hit a 0-sheet
            // attach it to the one sheet and reciprocally

            // attach the 0-sheet to the 1-sheet
            originalData_.sheet1List_[sheet1Id].sheet0List_.push_back(
              originalData_.vertex2sheet0_[vertexId0]);

            // attach the 1-sheet to the 0-sheet
            originalData_.sheet0List_[originalData_.vertex2sheet0_[vertexId0]]
              .sheet1List_.push_back(sheet1Id);
          }

          if(originalData_.vertex2sheet0_[vertexId1]
             >= static_cast<SimplexId>(originalData_.sheet0List_.size())) {
            // WEIRD BUG after multiple runs
            originalData_.vertex2sheet0_[vertexId1] = -1;
          }
          // now do the same for the other vertex
          // make sure these are not already visited 0-sheets
          if(originalData_.vertex2sheet0_[vertexId1] == -1) {

            // inspect neighbor edges
            SimplexId neighborNumber = 0;
            SimplexId nextEdgeId = -1;

            SimplexId vertexEdgeNumber
              = triangulation.getVertexEdgeNumber(vertexId1);

            for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

              SimplexId vertexEdgeId = -1;
              triangulation.getVertexEdge(vertexId1, j, vertexEdgeId);

              if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
                neighborNumber++;
                if(vertexEdgeId != edgeId
                   /*&&(saddleEdge[(*vertexEdgeList_)[vertexId1][j]] == 1)*/) {
                  nextEdgeId = vertexEdgeId;
                }
              }
            }

            if((neighborNumber == 2) && (nextEdgeId != -1)) {
              // this is a regular vertex and we can add the next edge to the
              // queue if it's not been visited before.
              if(!visitedEdges[nextEdgeId])
                edgeQueue.push(nextEdgeId);
            } else {
              // we hit a 0-sheet that hasn't been visited before.
              // let's create it

              // mark up the segmentation
              originalData_.vertex2sheet0_[vertexId1]
                = originalData_.sheet0List_.size();

              originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                               + 1);
              originalData_.sheet0List_.back().vertexId_ = vertexId1;
              originalData_.sheet0List_.back().type_ = 1;
              originalData_.sheet0List_.back().pruned_ = false;
            }
          }

          if(originalData_.vertex2sheet0_[vertexId1] != -1) {
            // we hit a 0-sheet
            // attach it to the one sheet and reciprocally

            // attach the 0-sheet to the 1-sheet
            originalData_.sheet1List_[sheet1Id].sheet0List_.push_back(
              originalData_.vertex2sheet0_[vertexId1]);

            // attach the 1-sheet to the 0-sheet
            originalData_.sheet0List_[originalData_.vertex2sheet0_[vertexId1]]
              .sheet1List_.push_back(sheet1Id);
          }

          visitedEdges[edgeId] = true;
        }

      } while(edgeQueue.size());
    }
  }

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#1-sheets", std::to_string(originalData_.sheet1List_.size())},
    {"#0-sheets", std::to_string(originalData_.sheet0List_.size())}});
  this->printMsg(
    "Extracted 0- and 1-sheets", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::compute3sheet(
  const SimplexId &vertexId,
  const std::vector<std::vector<std::vector<SimplexId>>> &tetTriangles,
  const triangulationType &triangulation) {

  SimplexId sheetId = originalData_.sheet3List_.size();
  originalData_.sheet3List_.resize(originalData_.sheet3List_.size() + 1);
  originalData_.sheet3List_.back().pruned_ = false;
  originalData_.sheet3List_.back().preMerger_ = -1;
  originalData_.sheet3List_.back().Id_ = sheetId;

  std::queue<SimplexId> vertexQueue;
  vertexQueue.push(vertexId);

  do {

    SimplexId localVertexId = vertexQueue.front();
    vertexQueue.pop();

    if(originalData_.vertex2sheet3_[localVertexId] == -1) {
      // not visited yet

      originalData_.sheet3List_.back().vertexList_.push_back(localVertexId);
      originalData_.vertex2sheet3_[localVertexId] = sheetId;

      SimplexId vertexStarNumber
        = triangulation.getVertexStarNumber(localVertexId);

      for(SimplexId i = 0; i < vertexStarNumber; i++) {
        SimplexId tetId = -1;
        triangulation.getVertexStar(localVertexId, i, tetId);

        if(tetTriangles[tetId].empty()) {
          //             if(originalData_.tet2sheet3_[tetId] == -1){
          //               originalData_.tet2sheet3_[tetId] = sheetId;
          //               originalData_.sheet3List_.back().tetList_.push_back(tetId);
          //             }
          for(int j = 0; j < 4; j++) {

            SimplexId tetVertexId = -1;
            triangulation.getCellVertex(tetId, j, tetVertexId);

            if(originalData_.vertex2sheet3_[tetVertexId] == -1) {
              vertexQueue.push(tetVertexId);
            }
          }
        } else {
          // fiber surface in there
          for(int j = 0; j < 4; j++) {
            SimplexId otherVertexId = -1;
            triangulation.getCellVertex(tetId, j, otherVertexId);
            if((otherVertexId != localVertexId)
               && (originalData_.vertex2sheet3_[otherVertexId] == -1)) {
              // we need to see if the edge <localVertexId, otherVertexId> is
              // cut by a fiber surface triangle or not.

              bool isCut = false;
              for(size_t k = 0; k < tetTriangles[tetId].size(); k++) {
                SimplexId l = 0, m = 0, n = 0;
                l = tetTriangles[tetId][k][0];
                m = tetTriangles[tetId][k][1];
                n = tetTriangles[tetId][k][2];

                for(int p = 0; p < 3; p++) {
                  std::pair<SimplexId, SimplexId> meshEdge;

                  if(fiberSurfaceVertexList_.size()) {
                    // the fiber surfaces have been merged
                    meshEdge
                      = fiberSurfaceVertexList_[originalData_.sheet2List_[l]
                                                  .triangleList_[m][n]
                                                  .vertexIds_[p]]
                          .meshEdge_;
                  } else {
                    // the fiber surfaces have not been merged
                    meshEdge = originalData_.sheet2List_[l]
                                 .vertexList_[m][originalData_.sheet2List_[l]
                                                   .triangleList_[m][n]
                                                   .vertexIds_[p]]
                                 .meshEdge_;
                  }

                  if(((meshEdge.first == localVertexId)
                      && (meshEdge.second == otherVertexId))
                     || ((meshEdge.second == localVertexId)
                         && (meshEdge.first == otherVertexId))) {
                    isCut = true;
                    break;
                  }
                }

                if(isCut)
                  break;
              }

              if(!isCut) {
                // add the vertex to the queue
                vertexQueue.push(otherVertexId);
              }
            }
          }
        }
      }
    }

  } while(vertexQueue.size());

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::compute3sheets(
  std::vector<std::vector<std::vector<SimplexId>>> &tetTriangles,
  const triangulationType &triangulation) {

  Timer t;

  tetTriangles.resize(tetNumber_);

  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {
    for(size_t j = 0; j < originalData_.sheet2List_[i].triangleList_.size();
        j++) {
      for(size_t k = 0;
          k < originalData_.sheet2List_[i].triangleList_[j].size(); k++) {

        SimplexId tetId
          = originalData_.sheet2List_[i].triangleList_[j][k].tetId_;

        std::vector<SimplexId> triangle(3);
        triangle[0] = i;
        triangle[1] = j;
        triangle[2] = k;

        tetTriangles[tetId].push_back(triangle);
      }
    }
  }

  // mark all the jacobi edge vertices
  for(size_t i = 0; i < originalData_.sheet1List_.size(); i++) {
    for(size_t j = 0; j < originalData_.sheet1List_[i].edgeList_.size(); j++) {

      SimplexId edgeId = originalData_.sheet1List_[i].edgeList_[j];

      SimplexId vertexId0 = -1, vertexId1 = -1;
      triangulation.getEdgeVertex(edgeId, 0, vertexId0);
      triangulation.getEdgeVertex(edgeId, 1, vertexId1);

      originalData_.vertex2sheet3_[vertexId0] = -2 - i;
      originalData_.vertex2sheet3_[vertexId1] = -2 - i;
    }
  }

  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(originalData_.vertex2sheet3_[i] == -1) {
      compute3sheet(i, tetTriangles, triangulation);
    }
  }

  // for 3-sheet expansion
  std::vector<std::vector<std::pair<SimplexId, bool>>> neighborList(
    originalData_.sheet3List_.size());
  // end of 3-sheet expansion

  // add the tets in parallel
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
    for(size_t j = 0; j < originalData_.sheet3List_[i].vertexList_.size();
        j++) {

      SimplexId vertexId = originalData_.sheet3List_[i].vertexList_[j];
      SimplexId sheetId = originalData_.vertex2sheet3_[vertexId];

      SimplexId vertexStarNumber = triangulation.getVertexStarNumber(vertexId);

      for(SimplexId k = 0; k < vertexStarNumber; k++) {
        SimplexId tetId = -1;
        triangulation.getVertexStar(vertexId, k, tetId);
        if(tetTriangles[tetId].empty()) {
          if(originalData_.tet2sheet3_[tetId] == -1) {
            originalData_.tet2sheet3_[tetId] = i;
            originalData_.sheet3List_[i].tetList_.push_back(tetId);
          }
        } else {

          // expending here 3-sheets.
          if(expand3sheets_) {
            for(int l = 0; l < 4; l++) {
              SimplexId otherVertexId = -1;

              triangulation.getCellVertex(tetId, l, otherVertexId);

              if(vertexId != otherVertexId) {

                SimplexId otherSheetId
                  = originalData_.vertex2sheet3_[otherVertexId];

                if((sheetId != otherSheetId) && (otherSheetId >= 0)) {

                  bool inThere = false;
                  for(size_t m = 0; m < neighborList[sheetId].size(); m++) {
                    if(neighborList[sheetId][m].first == otherSheetId) {
                      inThere = true;
                      break;
                    }
                  }

                  if(!inThere) {
                    neighborList[sheetId].push_back(
                      std::pair<SimplexId, bool>(otherSheetId, true));
                  }

                  for(size_t m = 0; m < tetTriangles[tetId].size(); m++) {

                    // see if this guy is a saddle
                    SimplexId x = tetTriangles[tetId][m][0];
                    SimplexId y = tetTriangles[tetId][m][1];
                    SimplexId z = tetTriangles[tetId][m][2];

                    FiberSurface::Triangle *tr
                      = &(originalData_.sheet2List_[x].triangleList_[y][z]);

                    bool cuttingTriangle = false;
                    for(int n = 0; n < 3; n++) {
                      FiberSurface::Vertex *v
                        = &(originalData_.sheet2List_[x]
                              .vertexList_[y][tr->vertexIds_[n]]);

                      if(((v->meshEdge_.first == vertexId)
                          && (v->meshEdge_.second == otherVertexId))
                         || ((v->meshEdge_.second == vertexId)
                             && (v->meshEdge_.first == otherVertexId))) {

                        cuttingTriangle = true;
                        break;
                      }
                    }

                    if(cuttingTriangle) {

                      SimplexId polygonId = tr->polygonEdgeId_;
                      SimplexId edgeId = jacobi2edges_[polygonId];
                      if(originalData_.edgeTypes_[edgeId] == 1) {
                        // this is a saddle Jacobi edge

                        inThere = false;
                        for(size_t n = 0; n < neighborList[sheetId].size();
                            n++) {

                          if(neighborList[sheetId][n].first == otherSheetId) {
                            if(neighborList[sheetId][n].second == true) {
                              neighborList[sheetId][n].second = false;
                            }
                            inThere = true;
                            break;
                          }
                        }

                        if(!inThere) {
                          neighborList[sheetId].push_back(
                            std::pair<SimplexId, bool>(otherSheetId, false));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          // end of the 3-sheet expansion
        }
      }
    }
  }

  SimplexId totalSheetNumber = 0;

  if(expand3sheets_) {

    totalSheetNumber = originalData_.sheet3List_.size();

    // expending sheets
    for(size_t i = 0; i < originalData_.sheet3List_.size(); i++) {
      if(originalData_.sheet3List_[i].pruned_ == false) {

        for(size_t j = 0; j < neighborList[i].size(); j++) {

          if(neighborList[i][j].second) {

            bool isForbidden = false;
            SimplexId neighborId = neighborList[i][j].first;

            // get the sheet where this guys has been merged
            while(originalData_.sheet3List_[neighborId].preMerger_ != -1) {

              neighborId = originalData_.sheet3List_[neighborId].preMerger_;
            }

            // make sure that no forbidden neighbor has been merged in the
            // candidate neighbor
            for(size_t k = 0;
                k
                < originalData_.sheet3List_[neighborId].preMergedSheets_.size();
                k++) {

              SimplexId subNeighborId
                = originalData_.sheet3List_[neighborId].preMergedSheets_[k];

              for(size_t l = 0; l < neighborList[i].size(); l++) {
                if((neighborList[i][l].first == subNeighborId)
                   && (!neighborList[i][l].second)) {
                  isForbidden = true;
                  break;
                }
              }
              if(isForbidden)
                break;
            }

            // make sure that neighborId is not a candidate for a merge with a
            // sheet that is forbidden for i
            if(!isForbidden) {
              for(size_t k = 0; k < neighborList[i].size(); k++) {
                if(!neighborList[i][k].second) {
                  SimplexId forbiddenNeighbor = neighborList[i][k].first;

                  // make sure forbiddenNeighbor is not a valid merger for
                  // neighborId
                  for(size_t l = 0; l < neighborList[neighborId].size(); l++) {
                    if((forbiddenNeighbor == neighborList[neighborId][l].first)
                       && (neighborList[neighborId][l].second)) {
                      isForbidden = true;
                      break;
                    }
                  }
                  if(isForbidden)
                    break;
                }
              }
            }

            if((neighborId != static_cast<SimplexId>(i)) && (!isForbidden)
               && (originalData_.sheet3List_[neighborId].pruned_ == false)
               && (originalData_.sheet3List_[neighborId].vertexList_.size()
                   > originalData_.sheet3List_[i].vertexList_.size())) {

              // these can be merged
              preMergeSheets(i, neighborId);
              totalSheetNumber--;
              break;
            }
          }
        }
      }
    }
  } else {
    totalSheetNumber = originalData_.sheet3List_.size();
  }
  // end of the 3-sheet expansion

  this->printMsg("Computed " + std::to_string(totalSheetNumber) + " 3-sheets",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::connectSheets(const triangulationType &triangulation) {

  Timer t;

  for(size_t i = 0; i < originalData_.sheet2List_.size(); i++) {
    for(size_t j = 0; j < originalData_.sheet2List_[i].triangleList_.size();
        j++) {

      for(size_t k = 0;
          k < originalData_.sheet2List_[i].triangleList_[j].size(); k++) {

        SimplexId tetId
          = originalData_.sheet2List_[i].triangleList_[j][k].tetId_;

        for(int l = 0; l < 4; l++) {
          SimplexId vertexId = -1;
          triangulation.getCellVertex(tetId, l, vertexId);

          SimplexId sheet3Id = originalData_.vertex2sheet3_[vertexId];

          if(sheet3Id >= 0) {
            connect3sheetTo2sheet(originalData_, sheet3Id, i);
          }
        }
      }
    }
  }

  // connect 3-sheets together
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(originalData_.vertex2sheet3_[i] >= 0) {

      SimplexId vertexEdgeNumber = triangulation.getVertexEdgeNumber(i);

      for(SimplexId j = 0; j < vertexEdgeNumber; j++) {
        SimplexId edgeId = -1;
        triangulation.getVertexEdge(i, j, edgeId);
        SimplexId otherVertexId = -1;
        triangulation.getEdgeVertex(edgeId, 0, otherVertexId);
        if(otherVertexId == i) {
          triangulation.getEdgeVertex(edgeId, 1, otherVertexId);
        }

        if(originalData_.vertex2sheet3_[otherVertexId] >= 0) {
          if(originalData_.vertex2sheet3_[otherVertexId]
             != originalData_.vertex2sheet3_[i]) {

            connect3sheetTo3sheet(originalData_,
                                  originalData_.vertex2sheet3_[i],
                                  originalData_.vertex2sheet3_[otherVertexId]);
          }
        }

        if(originalData_.vertex2sheet0_[otherVertexId] != -1) {
          connect3sheetTo0sheet(originalData_, originalData_.vertex2sheet3_[i],
                                originalData_.vertex2sheet0_[otherVertexId]);
        }

        if(originalData_.vertex2sheet3_[otherVertexId] < -1) {
          SimplexId sheet1Id = -2 - originalData_.vertex2sheet3_[otherVertexId];
          connect3sheetTo1sheet(
            originalData_, originalData_.vertex2sheet3_[i], sheet1Id);
        }
      }
    }
  }

  this->printMsg("Sheet connectivity established.");

  printConnectivity(originalData_);

  hasConnectedSheets_ = true;

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::disconnect3sheetFrom1sheet(
  ReebSpaceData &data,
  const SimplexId &sheet3Id,
  const SimplexId &sheet1Id,
  const SimplexId &biggerId,
  const triangulationType &triangulation) {

  std::vector<SimplexId> newList;

  newList.reserve(data.sheet1List_[sheet1Id].sheet3List_.size());
  for(size_t i = 0; i < data.sheet1List_[sheet1Id].sheet3List_.size(); i++) {
    SimplexId other3SheetId = data.sheet1List_[sheet1Id].sheet3List_[i];
    if((other3SheetId != sheet3Id) && (!data.sheet3List_[other3SheetId].pruned_)
       && (data.sheet3List_[other3SheetId].tetList_.size())) {
      newList.push_back(data.sheet1List_[sheet1Id].sheet3List_[i]);
    }
  }

  data.sheet1List_[sheet1Id].sheet3List_ = newList;

  if((data.sheet1List_[sheet1Id].hasSaddleEdges_)
     && (data.sheet1List_[sheet1Id].sheet3List_.size() == 1)) {
    // this guy is no longer separting any body.
    data.sheet1List_[sheet1Id].pruned_ = true;
    data.sheet2List_[sheet1Id].pruned_ = true;

    // update the segmentation
    for(size_t i = 0; i < data.sheet1List_[sheet1Id].edgeList_.size(); i++) {

      SimplexId vertexId = -1;
      triangulation.getEdgeVertex(
        data.sheet1List_[sheet1Id].edgeList_[i], 0, vertexId);
      data.vertex2sheet3_[vertexId] = biggerId;

      vertexId = -1;
      triangulation.getEdgeVertex(
        data.sheet1List_[sheet1Id].edgeList_[i], 1, vertexId);
      data.vertex2sheet3_[vertexId] = biggerId;
    }

    for(size_t i = 0; i < data.sheet1List_[sheet1Id].sheet0List_.size(); i++) {

      disconnect1sheetFrom0sheet(
        data, sheet1Id, data.sheet1List_[sheet1Id].sheet0List_[i], biggerId);
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::flush(const triangulationType &triangulation) {

  totalArea_ = -1;
  totalVolume_ = -1;
  totalHyperVolume_ = -1;
  hasConnectedSheets_ = false;

  // store the segmentation for later purpose
  originalData_.vertex2sheet0_.resize(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++)
    originalData_.vertex2sheet0_[i] = -1;

  originalData_.vertex2sheet3_.resize(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++)
    originalData_.vertex2sheet3_[i] = -1;

  originalData_.edge2sheet1_.resize(triangulation.getNumberOfEdges(), -1);
  originalData_.edgeTypes_.resize(triangulation.getNumberOfEdges(), -1);

  originalData_.tet2sheet3_.resize(tetNumber_, -1);

  jacobi2edges_.clear();
  jacobiSetEdges_.clear();

  originalData_.sheet0List_.clear();
  originalData_.sheet1List_.clear();
  originalData_.sheet2List_.clear();
  originalData_.sheet3List_.clear();

  fiberSurfaceVertexList_.clear();

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::mergeSheets(const SimplexId &smallerId,
                                const SimplexId &biggerId,
                                const triangulationType &triangulation) {

  // 1. add the vertices and tets of smaller to bigger
  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].vertexList_.size();
      i++) {
    SimplexId vertexId = currentData_.sheet3List_[smallerId].vertexList_[i];
    currentData_.sheet3List_[biggerId].vertexList_.push_back(vertexId);
    currentData_.vertex2sheet3_[vertexId] = biggerId;
  }
  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].tetList_.size();
      i++) {
    SimplexId tetId = currentData_.sheet3List_[smallerId].tetList_[i];
    currentData_.sheet3List_[biggerId].tetList_.push_back(tetId);
    currentData_.tet2sheet3_[tetId] = biggerId;
  }

  // 2. update bigger's score and re-insert it in the candidate list
  currentData_.sheet3List_[biggerId].domainVolume_
    += currentData_.sheet3List_[smallerId].domainVolume_;
  currentData_.sheet3List_[biggerId].rangeArea_
    += currentData_.sheet3List_[smallerId].rangeArea_;
  currentData_.sheet3List_[biggerId].hyperVolume_
    += currentData_.sheet3List_[smallerId].hyperVolume_;

  // 3. add smaller's connections to bigger's
  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet3List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet3List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo3sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet2List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet2List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo2sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet1List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet1List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo1sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet0List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet0List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo0sheet(currentData_, biggerId, otherSheetId);
    }
  }

  currentData_.sheet3List_[smallerId].pruned_ = true;

  // 4. disconnect smaller from all of its connections.
  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet3List_.size();
      i++) {
    disconnect3sheetFrom3sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet3List_[i]);
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet2List_.size();
      i++) {
    disconnect3sheetFrom2sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet2List_[i]);
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet1List_.size();
      i++) {
    disconnect3sheetFrom1sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet1List_[i], biggerId,
      triangulation);
  }

  for(size_t i = 0; i < currentData_.sheet3List_[smallerId].sheet0List_.size();
      i++) {
    disconnect3sheetFrom0sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet0List_[i]);
  }

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::simplifySheet(const SimplexId &sheetId,
                                  const SimplificationCriterion &criterion,
                                  const triangulationType &triangulation) {

  SimplexId candidateId = -1;
  double maximumScore = -1;

  // see the adjacent 3-sheets
  for(size_t i = 0; i < currentData_.sheet3List_[sheetId].sheet3List_.size();
      i++) {
    SimplexId otherSheetId = currentData_.sheet3List_[sheetId].sheet3List_[i];
    if((!currentData_.sheet3List_[otherSheetId].pruned_)
       && (sheetId != otherSheetId)) {

      double otherScore = 0;

      switch(criterion) {
        case SimplificationCriterion::domainVolume:
          otherScore = currentData_.sheet3List_[otherSheetId].domainVolume_;
          break;
        case SimplificationCriterion::rangeArea:
          otherScore = currentData_.sheet3List_[otherSheetId].rangeArea_;
          break;
        case SimplificationCriterion::hyperVolume:
          otherScore = currentData_.sheet3List_[otherSheetId].hyperVolume_;
          break;
      }

      if((maximumScore < 0) || (otherScore > maximumScore)) {
        candidateId = otherSheetId;
        maximumScore = otherScore;
      }
    }
  }

  // see adjacent 3-sheets along 1-sheets
  //   for(int i = 0;
  //     i < (int) currentData_.sheet3List_[sheetId].sheet1List_.size(); i++){
  //
  //     int sheet1Id = currentData_.sheet3List_[sheetId].sheet1List_[i];
  //
  //     for(int j = 0;
  //       j < (int) currentData_.sheet1List_[sheet1Id].sheet3List_.size();
  //       j++){
  //
  //       int otherSheetId = currentData_.sheet1List_[sheet1Id].sheet3List_[j];
  //       if((otherSheetId != sheetId)
  //         &&(!currentData_.sheet3List_[otherSheetId].pruned_)){
  //
  //         double otherScore = 0;
  //
  //         switch(criterion){
  //           case domainVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].domainVolume_;
  //             break;
  //           case rangeArea:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].rangeArea_;
  //             break;
  //           case hyperVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].hyperVolume_;
  //             break;
  //         }
  //
  //         if((maximumScore < 0)||(otherScore > maximumScore)){
  //           candidateId = otherSheetId;
  //           maximumScore = otherScore;
  //         }
  //       }
  //     }
  //   }

  // see adjacent 3-sheets on 0-sheets
  //   for(int i = 0;
  //     i < (int) currentData_.sheet3List_[sheetId].sheet0List_.size(); i++){
  //
  //     int sheet0Id = currentData_.sheet3List_[sheetId].sheet0List_[i];
  //
  //     for(int j = 0;
  //       j < (int) currentData_.sheet0List_[sheet0Id].sheet3List_.size();
  //       j++){
  //
  //       int otherSheetId = currentData_.sheet0List_[sheet0Id].sheet3List_[j];
  //       if((otherSheetId != sheetId)
  //         &&(!currentData_.sheet3List_[otherSheetId].pruned_)){
  //
  //         double otherScore = 0;
  //
  //         switch(criterion){
  //           case domainVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].domainVolume_;
  //             break;
  //           case rangeArea:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].rangeArea_;
  //             break;
  //           case hyperVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].hyperVolume_;
  //             break;
  //         }
  //
  //         if((maximumScore < 0)||(otherScore > maximumScore)){
  //           candidateId = otherSheetId;
  //           maximumScore = otherScore;
  //         }
  //       }
  //     }
  //   }

  // mark as pruned if so
  if(candidateId != -1) {
    mergeSheets(sheetId, candidateId, triangulation);
  }

  currentData_.sheet3List_[sheetId].pruned_ = true;

  return 0;
}

template <typename triangulationType>
int ttk::ReebSpace::simplifySheets(
  const double &simplificationThreshold,
  const SimplificationCriterion &simplificationCriterion,
  const triangulationType &triangulation) {

  Timer t;

  if(!currentData_.sheet3List_.size())
    return -1;

  SimplexId simplifiedSheets = 0;
  double lastThreshold = -1;

  for(size_t it = 0; it < originalData_.sheet3List_.size(); it++) {

    // do while, avoiding infinite loop

    double minValue = -1;
    SimplexId minId = -1;

    for(size_t i = 0; i < currentData_.sheet3List_.size(); i++) {
      if(!currentData_.sheet3List_[i].pruned_) {

        double value = 0;
        switch(simplificationCriterion) {

          case SimplificationCriterion::domainVolume:
            value = currentData_.sheet3List_[i].domainVolume_ / totalVolume_;
            break;

          case SimplificationCriterion::rangeArea:
            value = currentData_.sheet3List_[i].rangeArea_ / totalArea_;
            break;

          case SimplificationCriterion::hyperVolume:
            value
              = currentData_.sheet3List_[i].hyperVolume_ / totalHyperVolume_;
            break;
        }

        if((minId == -1) || (value < minValue)) {
          minValue = value;
          minId = i;
        }
      }
    }

    if((minId != -1) && (minValue < simplificationThreshold)) {
      simplifySheet(minId, simplificationCriterion, triangulation);
      simplifiedSheets++;
      lastThreshold = minValue;
    } else {
      break;
    }
  }

  currentData_.simplificationThreshold_ = simplificationThreshold;
  currentData_.simplificationCriterion_ = simplificationCriterion;

  SimplexId simplificationId = 0;
  for(size_t i = 0; i < currentData_.sheet3List_.size(); i++) {
    if(!currentData_.sheet3List_[i].pruned_) {
      currentData_.sheet3List_[i].simplificationId_ = simplificationId;
      simplificationId++;
    }
  }
  for(size_t i = 0; i < currentData_.sheet3List_.size(); i++) {
    if(currentData_.sheet3List_[i].pruned_) {
      // find where it merged
      if(currentData_.sheet3List_[i].vertexList_.size()) {
        SimplexId vertexId = currentData_.sheet3List_[i].vertexList_[0];
        SimplexId sheetId = currentData_.vertex2sheet3_[vertexId];
        if(sheetId != static_cast<SimplexId>(i)) {
          currentData_.sheet3List_[i].simplificationId_
            = currentData_.sheet3List_[sheetId].simplificationId_;
        }
      }
    }
  }

  // orphans
  for(size_t i = 0; i < currentData_.sheet1List_.size(); i++) {
    if((!currentData_.sheet1List_[i].pruned_)
       && (currentData_.sheet1List_[i].hasSaddleEdges_)) {
      SimplexId nonSimplified = 0;
      for(size_t j = 0; j < currentData_.sheet1List_[i].sheet3List_.size();
          j++) {
        SimplexId sheet3Id = currentData_.sheet1List_[i].sheet3List_[j];

        if(!currentData_.sheet3List_[sheet3Id].pruned_) {
          nonSimplified++;
        }
      }

      if(nonSimplified < 2) {
        currentData_.sheet1List_[i].pruned_ = true;
      }
    }
  }

  // TODO: update segmentation for 1-sheets and 0-sheets?...

  printConnectivity(currentData_);

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#3-sheets simplified", std::to_string(simplifiedSheets)},
    {"Last 3-sheet threshold", std::to_string(lastThreshold)},
    {"#3-sheets left", std::to_string(simplificationId)}});

  this->printMsg(
    "3-sheets simplified", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
