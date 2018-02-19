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

#ifndef _REEBSPACE_H
#define _REEBSPACE_H

// standard includes

// base code includes
#include                  <FiberSurface.h>
#include                  <Geometry.h>
#include                  <JacobiSet.h>
#include                  <Triangulation.h>
// to add FiberSurface.h
#include                  <Wrapper.h>

#include                  <map>
#include                  <set>

namespace ttk{
  
  class ReebSpace : public Debug{

    public:
        
      enum SimplificationCriterion{
        domainVolume, // 0
        rangeArea, // 1
        hyperVolume //2
      };
      
      class Sheet0{
        
        public: 
          
          int             vertexId_, pruned_;
          // 0: extrema-sheet, 1: saddle-sheet
          char            type_;
          vector<int>     sheet1List_;
          vector<int>     sheet3List_;
      };
      
      class Sheet1{
        
        public:
       
          bool            hasSaddleEdges_, pruned_;
          vector<int>     edgeList_;
          // 0: extrema-sheet, 1: saddle-sheet (can be inferred by the number
          // of 2-sheets attached to it?)
          vector<int>     sheet0List_;
          // NB: the corresponding 2-sheet should have the same global 
          // indentifier.
          vector<int>     sheet3List_;
      };
      
      class Sheet2{
        
        public:
        
          bool            pruned_;
          int             sheet1Id_;
          // for each point, x, y, z coordinates
          // this is how coincident point will be merged afterwards
          // this list is meant to be temporary. after the merge, we'll use a 
          // global list
          vector<vector<FiberSurface::Vertex> >
                          vertexList_;
          vector<vector<FiberSurface::Triangle> >
                          triangleList_;
          vector<int>     sheet3List_;
      };
      
      class Sheet3{
        
        public:
        
          int             Id_, simplificationId_, preMerger_;
          bool            pruned_;
          double          domainVolume_, rangeArea_, hyperVolume_;
          vector<int>     vertexList_;
          vector<int>     tetList_;
          vector<int>     sheet0List_;
          vector<int>     sheet1List_;
          vector<int>     sheet2List_;
          vector<int>     sheet3List_;
          vector<int>     preMergedSheets_;
      };
      
      ReebSpace();
      
      ~ReebSpace();
      
      inline bool empty() const{
        return (currentData_.vertex2sheet0_.size() == 0);
      }
          
      template <class dataTypeU, class dataTypeV> inline int execute();
      
      inline const Sheet0* get0sheet(const int &sheetId) const{
        
#ifndef TTK_ENABLE_KAMIKAZE
        if((sheetId < 0)||(sheetId >= (int) currentData_.sheet0List_.size()))
          return NULL;
#endif
        
        return &(currentData_.sheet0List_[sheetId]);
      }
      
      inline const Sheet1* get1sheet(const int &sheetId) const{
#ifndef TTK_ENABLE_KAMIKAZE
        if((sheetId < 0)||(sheetId >= (int) currentData_.sheet1List_.size()))
          return NULL;
#endif
        
        return &(currentData_.sheet1List_[sheetId]); 
      }
      
      // warning, these are the originals
      inline const Sheet2* get2sheet(const int &sheetId) const{
#ifndef TTK_ENABLE_KAMIKAZE
        if((sheetId < 0)||(sheetId >= (int) originalData_.sheet2List_.size()))
          return NULL;
#endif
        
        return &(originalData_.sheet2List_[sheetId]);
      }
      
      inline const Sheet3* get3sheet(const int &sheetId) const{
#ifndef TTK_ENABLE_KAMIKAZE
        if((sheetId < 0)||(sheetId >= (int) currentData_.sheet3List_.size()))
          return NULL;
#endif
        
        return &(currentData_.sheet3List_[sheetId]);
      }
      
      inline const vector<int>* get0sheetSegmentation() const{
        return &currentData_.vertex2sheet0_;
      }
      
      const vector<int>* get1sheetSegmentation() const{
        return &currentData_.edge2sheet1_;
      }
     
      const vector<int>* get3sheetVertexSegmentation() const{
        return &currentData_.vertex2sheet3_;
      }
      
      const vector<int>* get3sheetTetSegmentation() const{
        return &currentData_.tet2sheet3_;
      }
     
      const vector<int>* getEdgeTypes() const{
        return &currentData_.edgeTypes_;
      }
     
      inline const vector<FiberSurface::Vertex>* 
        getFiberSurfaceVertices() const{
        return &(fiberSurfaceVertexList_);
      }
      
      inline int getJacobi2Edge(const int &jacobiEdgeId) const{
        if((jacobiEdgeId < 0)
          ||(jacobiEdgeId >= (int) jacobi2edges_.size()))
          return -1;
        return jacobi2edges_[jacobiEdgeId];
      }
      
      inline int getNumberOf2sheets() const {
        return currentData_.sheet2List_.size();
      }
     
//       inline vector<long long int>* getSheetTriangulationCells(){
//         return &sheet3cells_;
//       }
//      
//       inline vector<float>* getSheetTriangulationPoints(){
//         return &sheet3points_;
//       }
     
      template <class dataTypeU, class dataTypeV>
        inline int perturbate(
          const dataTypeU &uEpsilon = pow(10, -DBL_DIG),
          const dataTypeV &vEpsilon = pow(10, -DBL_DIG)) const;
         
      inline int setExpand3Sheets(const bool &onOff){
        expand3sheets_ = onOff;
        return 0;
      }
      
      inline int setInputField(const void *uField, const void *vField){
        
        uField_ = uField;
        vField_ = vField;
        
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
        fiberSurface_.flushOctree();
#endif
        
        return 0;
      }
     
      inline bool setRangeDrivenOctree(const bool &onOff){
       
        if(onOff != withRangeDrivenOctree_){
          withRangeDrivenOctree_ = onOff;
          return true;
        }

        withRangeDrivenOctree_ = onOff;
        return false;
      }
      
      inline int setSosOffsetsU(vector<int> *sosOffsetsU){
        sosOffsetsU_ = sosOffsetsU;
        return 0;
      }
      
      inline int setSosOffsetsV(vector<int> *sosOffsetsV){
        sosOffsetsV_ = sosOffsetsV;
        return 0;
      }

      inline int setTetNumber(const int &tetNumber){
        tetNumber_ = tetNumber;
        return 0;
      }

      /// Set the number of vertices in the scalar field.
      /// \param vertexNumber Number of vertices in the data-set.
      /// \return Returns 0 upon success, negative values otherwise. 
      int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }
      
      // WARNING: if you plan to use the range driven octree, make sure
      // that you provided pointers to the u and v fields.
      template <class dataTypeU, class dataTypeV>
        inline int setupTriangulation(Triangulation *triangulation){
       
        triangulation_ = triangulation;
        
        if(triangulation_){
          
          triangulation_->preprocessVertexStars();
          triangulation_->preprocessEdges();
          triangulation_->preprocessVertexEdges();
          
          JacobiSet<dataTypeU, dataTypeV> jacobiSet;
          jacobiSet.setWrapper(wrapper_);
          
          // trigger the jacobiSet pre-processing on the triangulation.
          jacobiSet.setupTriangulation(triangulation_);
          
          // trigger the fiberSurface pre-processing on the triangulation.
          fiberSurface_.setWrapper(wrapper_);
          fiberSurface_.setInputField(uField_, vField_);
          fiberSurface_.setupTriangulation(triangulation_);
          // trigger the fiber surface precomputation
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
          if(withRangeDrivenOctree_)
            fiberSurface_.buildOctree<dataTypeU, dataTypeV>();
#endif
          
          vertexNumber_ = triangulation_->getNumberOfVertices();
          edgeNumber_ = triangulation_->getNumberOfEdges();
          tetNumber_ = triangulation_->getNumberOfCells();
        }
        
        return 0;
      }
      
      template <class dataTypeU, class dataTypeV>
        inline int simplify(const double &simplificationThreshold,
        const SimplificationCriterion &criterion = rangeArea);
      

    protected:
    
      class ReebSpaceData;
      
      int compute1sheetsOnly(const vector<pair<int, char> > &jacobiSet,
        vector<pair<int, int> > &jacobiSetClassification);
      
      int compute1sheets(const vector<pair<int, char> > &jacobiSet,
        vector<pair<int, int> > &jacobiSetClassification);
      
      template <class dataTypeU, class dataTypeV>
        inline int compute2sheets(
          const vector<pair<int, int> > &jacobiEdges);
        
      template <class dataTypeU, class dataTypeV>
        inline int compute2sheetChambers();
      
      int compute3sheet(const int &vertexId, 
        const vector<vector<vector<int> > > &tetTriangles);
      
      int compute3sheets(vector<vector<vector<int> > > &tetTriangles);
        
      template <class dataTypeU, class dataTypeV>
        inline int computeGeometricalMeasures(Sheet3 &sheet);
      
      int connect3sheetTo0sheet(
        ReebSpaceData &data,
        const int &sheet3Id, const int &sheet0Id);
       
      int connect3sheetTo1sheet(
        ReebSpaceData &data,
        const int &sheet3Id, const int &sheet1Id);
      
      int connect3sheetTo2sheet(
        ReebSpaceData &data,
        const int &sheet3Id, const int &sheet2Id);
      
      int connect3sheetTo3sheet(
        ReebSpaceData &data,
        const int &sheet3Id, const int &otherSheet3Id);
      
      int connectSheets();
        
      int disconnect1sheetFrom0sheet(ReebSpaceData &data,
        const int &sheet1Id, const int &sheet0Id, const int &biggerId);
      
      int disconnect3sheetFrom0sheet(ReebSpaceData &data,
        const int &sheet3Id, const int &sheet0Id);
      
      int disconnect3sheetFrom1sheet(ReebSpaceData &data,
        const int &sheet3Id, const int &sheet1Id, const int &biggerId);
      
      int disconnect3sheetFrom2sheet(ReebSpaceData &data,
        const int &sheet3Id, const int &sheet2Id);
      
      int disconnect3sheetFrom3sheet(ReebSpaceData &data,
        const int &sheet3Id, const int &other3SheetId);
      
      int flush();
      
      int mergeSheets(const int &smallerId, const int &biggerId);
      
      int preMergeSheets(const int &sheetId0, const int &sheetId1);
      
      int prepareSimplification();
     
      int printConnectivity(ostream &stream, const ReebSpaceData &data) const;
      
      int simplifySheets(const double &simplificationThreshold,
        const SimplificationCriterion &simplificationCriterion);
     
      int simplifySheet(const int &sheetId, 
        const SimplificationCriterion &simplificationCriterion);
      
//       int triangulateTetrahedron(const int &tetId,
//         const vector<vector<int> > &triangles,
//         vector<long long int> &outputTets);
//       
//       int triangulateThreeSheets();
      
      int                   vertexNumber_, edgeNumber_, tetNumber_;
      double                totalArea_, totalVolume_, totalHyperVolume_;
      
      const void            *uField_, *vField_;
      vector<int>           *sosOffsetsU_, *sosOffsetsV_;
      
      // output segmentation 
      class ReebSpaceData{
        
        public: 
          
          SimplificationCriterion
                              simplificationCriterion_;
          double              simplificationThreshold_;
         
          vector<int>           edge2sheet1_;
          vector<int>           edgeTypes_;
          vector<int>           tet2sheet3_;
          vector<int>           vertex2sheet0_;
          vector<int>           vertex2sheet3_;
          
          // structure
          vector<Sheet0>        sheet0List_;
          vector<Sheet1>        sheet1List_;
          vector<Sheet2>        sheet2List_;
          vector<Sheet3>        sheet3List_;
          
//         vector<float>         sheet3points_;
//         vector<long long int> sheet3cells_;
      };
      
      bool                  hasConnectedSheets_, 
                            expand3sheets_,
                            withRangeDrivenOctree_;
      ReebSpaceData         originalData_, currentData_; 
      
      // information that does not get simplified
      vector<pair<int, char> >
                            jacobiSetEdges_;
      vector<int>           jacobi2edges_;
      
      FiberSurface          fiberSurface_;
      vector<FiberSurface::Vertex>
                            fiberSurfaceVertexList_;
                            
      Triangulation         *triangulation_;
  };
}

// if the package is not a template, comment the following line
// #include                  <ReebSpace.cpp>

template <class dataTypeU, class dataTypeV> 
  inline int ReebSpace::execute(){

  flush();
  
  Timer t;
  
  // 1) compute the jacobi set
  JacobiSet<dataTypeU, dataTypeV> jacobiSet;
  
  jacobiSet.setWrapper(wrapper_);
  jacobiSet.setupTriangulation(triangulation_);
  jacobiSet.setInputField(uField_, vField_);
  jacobiSet.setSosOffsetsU(sosOffsetsU_);
  jacobiSet.setSosOffsetsV(sosOffsetsV_);
  jacobiSet.execute(jacobiSetEdges_);
  
  // 2) compute the list saddle 1-sheets
  // + list of saddle 0-sheets
  vector<pair<int, int> > jacobiSetClassification;
  compute1sheetsOnly(
    jacobiSetEdges_, jacobiSetClassification);
  // at this stage, jacobiSetClassification contains the list of saddle edges 
  // along with their 1-sheet Id.
  
  compute2sheets<dataTypeU, dataTypeV>(jacobiSetClassification);
//   compute2sheetChambers<dataTypeU, dataTypeV>();
  
  vector<vector<vector<int> > > tetTriangles;
  compute3sheets(tetTriangles);
  
  {
    stringstream msg;
    msg << "[ReebSpace] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. TOTAL (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // post-process for further interaction
  if((totalArea_ == -1)||(totalVolume_ == -1)||(totalHyperVolume_ == -1)){
    
    Timer t;
    
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) originalData_.sheet3List_.size(); i++){
      computeGeometricalMeasures<dataTypeU, dataTypeV>(
        originalData_.sheet3List_[i]);
    }
    
    for(int i = 0; i < (int) originalData_.sheet3List_.size(); i++){
      totalArea_ += originalData_.sheet3List_[i].rangeArea_;
      totalVolume_ += originalData_.sheet3List_[i].domainVolume_;
      totalHyperVolume_ += originalData_.sheet3List_[i].hyperVolume_;
    }
    
    {
      stringstream msg;
      msg << "[ReebSpace] Geometrical measures computed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))" << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  
  fiberSurface_.finalize<dataTypeU, dataTypeV>();
  
  prepareSimplification();
  
  return 0;
}

template <class dataTypeU, class dataTypeV>
  inline int ReebSpace::compute2sheets(
    const vector<pair<int, int> > &jacobiEdges){

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
#endif
 
  Timer t;
  
  // at this point, they have exactly the same size
  originalData_.sheet2List_.resize(originalData_.sheet1List_.size());
  for(int i = 0; i < (int) originalData_.sheet2List_.size(); i++){
    originalData_.sheet2List_[i].sheet1Id_ = i;
    originalData_.sheet2List_[i].pruned_ = false;
    originalData_.sheet2List_[i].vertexList_.resize(
      originalData_.sheet1List_[
        originalData_.sheet2List_[i].sheet1Id_].edgeList_.size());
    originalData_.sheet2List_[i].triangleList_.resize(
      originalData_.sheet1List_[
        originalData_.sheet2List_[i].sheet1Id_].edgeList_.size());
    for(int j = 0; 
      j < (int) originalData_.sheet2List_[i].vertexList_.size(); j++){
      originalData_.sheet2List_[i].vertexList_[j].clear();
      originalData_.sheet2List_[i].triangleList_[j].clear();
    }
  }
  
  fiberSurface_.setWrapper(wrapper_);
  fiberSurface_.setGlobalVertexList(&fiberSurfaceVertexList_);
  fiberSurface_.setInputField(uField_, vField_);
  fiberSurface_.setupTriangulation(triangulation_);
  
  fiberSurface_.setPolygonEdgeNumber(jacobiEdges.size());
 
  vector<int> edge2polygonEdgeId(edgeNumber_, -1);
  jacobi2edges_.resize(jacobiEdges.size());
  
  for(int i = 0; i < (int) jacobiEdges.size(); i++){
    
    int edgeId = jacobiEdges[i].first;
    
    edge2polygonEdgeId[edgeId] = i;
    jacobi2edges_[i] = edgeId; 
  }
  
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) originalData_.sheet2List_.size(); i++){
    
    for(int j = 0; 
      j < (int) originalData_.sheet1List_[
        originalData_.sheet2List_[i].sheet1Id_].edgeList_.size(); j++){
     
      int edgeId = originalData_.sheet1List_[
        originalData_.sheet2List_[i].sheet1Id_].edgeList_[j];
        
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
  for(int i = 0; i < (int) jacobiEdges.size(); i++){
    
    int edgeId = jacobiEdges[i].first;
    
    pair<double, double> rangePoint0, rangePoint1;
    
    int vertexId0 = -1, vertexId1 = -1;
    triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
    triangulation_->getEdgeVertex(edgeId, 1, vertexId1);
    
    rangePoint0.first = ((dataTypeU *) uField_)[vertexId0];
    rangePoint0.second = ((dataTypeV *) vField_)[vertexId0];
      
    rangePoint1.first = ((dataTypeU *) uField_)[vertexId1];
    rangePoint1.second = ((dataTypeV *) vField_)[vertexId1];
      
    if(originalData_.edgeTypes_[edgeId] == 1){
      vector<int> edgeSeeds(triangulation_->getEdgeStarNumber(edgeId), -1);
      for(int j = 0; j < (int) edgeSeeds.size(); j++){
        triangulation_->getEdgeStar(edgeId, j, edgeSeeds[j]);
      }
      
      fiberSurface_.computeContour<dataTypeU, dataTypeV>(
        rangePoint0, rangePoint1,
        edgeSeeds,
        edge2polygonEdgeId[edgeId]);
    }
    else{
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
      if(withRangeDrivenOctree_){
        fiberSurface_.computeSurfaceWithOctree<dataTypeU, dataTypeV>(
          rangePoint0, rangePoint1, edge2polygonEdgeId[edgeId]);
      }
      else{
        fiberSurface_.computeSurface<dataTypeU, dataTypeV>(
          rangePoint0, rangePoint1, edge2polygonEdgeId[edgeId]);
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
//     vector<bool> inList(tetNumber_, false);
//     vector<int> seedList;
//     vector<pair<pair<double, double>, pair<double, double> > > edgeList;
//     vector<int> jacobiEdgeIdList;
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
  
  {
    stringstream msg;
    msg << "[ReebSpace] Fiber surfaces computed in "
      << t.getElapsedTime()
      << " s. overall (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <class dataTypeU, class dataTypeV>
  inline int ReebSpace::compute2sheetChambers(){

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
#endif
 
  {
    stringstream msg;
    msg << "[ReebSpace] Computing chambers' pre-images." << endl;
    msg << "[ReebSpace] This will take a LONG time." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  Timer t;
  
  // at this point, they have exactly the same size
  originalData_.sheet2List_.resize(threadNumber_);
  for(int i = 0; i < (int) originalData_.sheet2List_.size(); i++){
    originalData_.sheet2List_[i].sheet1Id_ = i;
    originalData_.sheet2List_[i].pruned_ = false;
    originalData_.sheet2List_[i].vertexList_.resize(1);
    originalData_.sheet2List_[i].triangleList_.resize(1);
  }
  
  fiberSurface_.setWrapper(wrapper_);
  fiberSurface_.setGlobalVertexList(&fiberSurfaceVertexList_);
  fiberSurface_.setInputField(uField_, vField_);
  fiberSurface_.setupTriangulation(triangulation_);
  fiberSurface_.setTetNumber(tetNumber_);
  fiberSurface_.setPolygonEdgeNumber(threadNumber_);
  
  for(int i = 0; i < (int) originalData_.sheet2List_.size(); i++){
    
    fiberSurface_.setTriangleList(i, 
      &(originalData_.sheet2List_[i].triangleList_[0]));
    fiberSurface_.setVertexList(i,
      &(originalData_.sheet2List_[i].vertexList_[0]));
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) edgeNumber_; i++){

 #ifdef TTK_ENABLE_OPENMP
    int threadId = omp_get_thread_num();
 #else
    int threadId = 0;
 #endif
    
    int edgeId = i;
    
    pair<double, double> rangePoint0, rangePoint1;
    
    int vertexId0 = -1, vertexId1 = -1;
    triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
    triangulation_->getEdgeVertex(edgeId, 1, vertexId1);
    
    rangePoint0.first = ((dataTypeU *) uField_)[vertexId0];
    rangePoint0.second = ((dataTypeV *) vField_)[vertexId0];
      
    rangePoint1.first = ((dataTypeU *) uField_)[vertexId1];
    rangePoint1.second = ((dataTypeV *) vField_)[vertexId1];
      
    fiberSurface_.computeSurface<dataTypeU, dataTypeV>(
      rangePoint0, rangePoint1, threadId);
    
    // clear the memory now.. otherwise we'll swap and never end
    originalData_.sheet2List_[threadId].triangleList_[0].clear();
    originalData_.sheet2List_[threadId].vertexList_[0].clear();
  }

  {
    stringstream msg;
    msg << "[ReebSpace] Chambers pre-image boundaries computed in "
      << t.getElapsedTime()
      << " s. overall (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV>
  inline int ReebSpace::computeGeometricalMeasures(Sheet3 &sheet){

  sheet.domainVolume_ = 0;
  sheet.rangeArea_ = 0;
  sheet.hyperVolume_ = 0;
  
  for(int i = 0; i < (int) sheet.tetList_.size(); i++){
    
    int tetId = sheet.tetList_[i];
    vector<pair<double, double> > domainBox, rangeBox;
    vector<vector<float> > domainPoints(4), rangePoints(4);
    
    for(int j = 0; j < 4; j++){
      domainPoints[j].resize(3);
      rangePoints[j].resize(2);
      
      int vertexId = -1;
      triangulation_->getCellVertex(tetId, j, vertexId);
      
      triangulation_->getVertexPoint(vertexId,
        domainPoints[j][0], domainPoints[j][1], domainPoints[j][2]);
      
      rangePoints[j][0] = ((dataTypeU *) uField_)[vertexId];
      rangePoints[j][1] = ((dataTypeV *) vField_)[vertexId];
    }
    
    Geometry::getBoundingBox(domainPoints, domainBox);
    Geometry::getBoundingBox(rangePoints, rangeBox);
    
    sheet.domainVolume_ += 
      (domainBox[0].second - domainBox[0].first)
      *(domainBox[1].second - domainBox[1].first)
      *(domainBox[2].second - domainBox[2].first);
    
    sheet.rangeArea_ += 
      (rangeBox[0].second - rangeBox[0].first)
      *(rangeBox[1].second - rangeBox[1].first);
    
  }
  
  if(sheet.domainVolume_){
    sheet.hyperVolume_ = sheet.rangeArea_/sheet.domainVolume_;
  }
  else{
    sheet.hyperVolume_ = 0;
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV>
  inline int ReebSpace::perturbate(
    const dataTypeU &uEpsilon, const dataTypeV &vEpsilon) const{

  JacobiSet<dataTypeU, dataTypeV> jacobiSet;
  jacobiSet.setWrapper(wrapper_);
  jacobiSet.setInputField(uField_, vField_);
  jacobiSet.setVertexNumber(vertexNumber_);
  jacobiSet.perturbate(uEpsilon, vEpsilon);
  
  return 0;
}

template <class dataTypeU, class dataTypeV>
  inline int ReebSpace::simplify(const double &simplificationThreshold, 
    const SimplificationCriterion &criterion){

  if((totalArea_ == -1)||(totalVolume_ == -1)||(totalHyperVolume_ == -1)){
    
    Timer t;
    
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) originalData_.sheet3List_.size(); i++){
      computeGeometricalMeasures<dataTypeU, dataTypeV>(
        originalData_.sheet3List_[i]);
    }
    
    for(int i = 0; i < (int) originalData_.sheet3List_.size(); i++){
      totalArea_ += originalData_.sheet3List_[i].rangeArea_;
      totalVolume_ += originalData_.sheet3List_[i].domainVolume_;
      totalHyperVolume_ += originalData_.sheet3List_[i].hyperVolume_;
    }
    
    {
      stringstream msg;
      msg << "[ReebSpace] Geometrical measures computed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))" << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  
  if(!hasConnectedSheets_){
    connectSheets();
    prepareSimplification();
  }
  
  {
    stringstream msg;
    msg << "[ReebSpace] Simplifying with criterion ";
    switch(criterion){
      case domainVolume:
        msg << "'Domain Volume'";
        break;
      case rangeArea:
        msg << "'Range Area'";
        break;
      case hyperVolume:
        msg << "'HyperVolume'";
        break;
    }
    msg << " at threshold " << simplificationThreshold << "." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  if(!((criterion == currentData_.simplificationCriterion_)
    &&(simplificationThreshold > currentData_.simplificationThreshold_))){
    prepareSimplification();
  }
  
  simplifySheets(simplificationThreshold, criterion);
  
  return 0;
}

#endif // REEBSPACE_H
