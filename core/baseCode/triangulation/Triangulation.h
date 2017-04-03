/// \ingroup baseCode
/// \class ttk::Triangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief Triangulation is a class that provides time and memory efficient 
/// traversal methods on triangulations of piecewise linear manifolds. It 
/// provides the following features:
///   -# Given a vertex, it provides: the list of edges that are connected to 
/// it, the list of its neighbors, its link, its star, etc.
///   -# Given an edge, it provides: its vertices, its star, etc.
///   -# Given a triangle, its provides: its vertices, its edges, etc.
///   -# Given a tetrahedron, its provides: its vertices, its edges, its 
/// neighbor tetrahedra, etc.
///   -# Given a triangulation, it provides: its list of vertices, edges, 
/// triangles and tetrahedra.
///
/// Triangulation supports both explicit and implicit triangulations: 
///   -# Explicit triangulations: Given a list of points and a list of cells,
/// Triangulation provides time efficient accesses (requiring adequate 
/// pre-processing, see the documentation furtherdown).
///   -# Implicit triangulations: Given a regular grid (origin, spacings and
/// dimensions), Triangulation will perform an implicit triangulation of the 
/// grid, enabling both time and memory efficient traversals of triangulations 
/// of regular grids.
///
/// Apart from pre-processes, Triangulation requires no memory overhead in
/// addition to the input data. 
/// 
/// \note
/// Only pre-process the information you need! See the documentation further 
/// down.
/// \sa vtkTriangulation

#ifndef _TRIANGULATION_H
#define _TRIANGULATION_H

// base code includes
#include                  <AbstractTriangulation.h>
#include                  <ImplicitTriangulation.h>
#include                  <ExplicitTriangulation.h>

namespace ttk{
  
  class Triangulation : public AbstractTriangulation{

    public:
        
      Triangulation();
      
      ~Triangulation();

      /// Reset the triangulation data-structures.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int clear(){
       
        if(abstractTriangulation_){
          return abstractTriangulation_->clear();
        }
        
        return 0;
      }
      
      /// Computes and displays the memory footprint of the data-structure.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int footprint() const {
        
        if(abstractTriangulation_){
          return abstractTriangulation_->footprint();
        }
        
        return 0;
      }
      
      /// Get the \p localEdgeId-th edge of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellEdges() needs to be called on this object prior to any 
      /// traversal, in a clearly distinct pre-processing step that involves no 
      /// traversal at all. An error will be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localEdgeId Input local edge identifier, 
      /// in [0, getCellEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellEdgeNumber()
      inline int getCellEdge(const int &cellId,
        const int &localEdgeId, int &edgeId) const{
         
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedCellEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellEdge query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellEdge(
          cellId, localEdgeId, edgeId);
      }

      /// Get the number of edges for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, preprocessCellEdges() 
      /// needs to be called on this object prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell edges.
      inline int getCellEdgeNumber(const int &cellId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        if(!abstractTriangulation_->hasPreprocessedCellEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellEdgeNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellEdgeNumber(cellId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of edges for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell edge list.
      inline const vector<vector<int> > *getCellEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedCellEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellEdges query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getCellEdges();
      }
      
      /// Get the \p localNeighborId-th cell neighbor of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localNeighborId Input local neighbor identifier, 
      /// in [0, getCellNeighborNumber()].
      /// \param neighborId Output global neighbor cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellNeighborNumber()
      inline int getCellNeighbor(const int &cellId,
        const int &localNeighborId, int &neighborId) const{
         
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedCellNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellNeighbor query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellNeighbor(
          cellId, localNeighborId, neighborId);
      }
      
      /// Get the number of cell neighbors for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell neighbors.
      inline int getCellNeighborNumber(const int &cellId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedCellNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellNeighborNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellNeighborNumber(cellId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of cell neighbors for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of neighbor cells for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell neighbor list.
      inline const vector<vector<int> > *getCellNeighbors(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedCellNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellNeighbors query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getCellNeighbors();
      }
      
      /// Get the \p localTriangleId-th triangle id of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getCellTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellTriangleNumber()
      inline int getCellTriangle(const int &cellId,
        const int &localTriangleId, int &triangleId) const{
         
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedCellTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellTriangle query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellTriangle(
          cellId, localTriangleId, triangleId);
      }
      
      /// Get the number of triangles for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell triangles.
      inline int getCellTriangleNumber(const int &cellId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedCellTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellTriangleNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getCellTriangleNumber(cellId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell triangle list.
      inline const vector<vector<int> > *getCellTriangles(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedCellTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "CellTriangles query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getCellTriangles();
      }
      
      /// Get the \p localVertexId-th vertex identifier of the \p cellId-th 
      /// cell.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \param cellId Input global cell identifier.
      /// \param localVertexId Input local vertex identifier,
      /// in [0, getCellVertexNumber()].
      /// \param vertexId Ouput global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellVertexNumber()
      inline int getCellVertex(const int &cellId,
        const int &localVertexId, int &vertexId) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
          
        return abstractTriangulation_->getCellVertex(
          cellId, localVertexId, vertexId);
      }
      
      /// Get the number of vertices in a cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \param cellId Input global cell identifier.
      /// \returns Number of vertices in the cell.
      inline int getCellVertexNumber(const int &cellId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        return abstractTriangulation_->getCellVertexNumber(cellId);
      }
      
      /// Get the dimensionality of the triangulation (this value is equal to 
      /// the dimension of the simplex with largest dimensionality).
      /// \return Returns the dimensionality of the triangulation.
      inline int getDimensionality() const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        return abstractTriangulation_->getDimensionality();
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges of the triangulation.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a pair of vertex identifiers.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge list.
      inline const vector<pair<int, int> > *getEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "Edges query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getEdges();
      }
      
      /// Get the \p localLinkId-th simplex of the link of the \p edgeId-th 
      /// edge.
      ///
      /// The output \p linkId refers in 2D to a vertex identifier and in 3D
      /// to an edge identifier.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getEdgeLinkNumber()].
      /// \param linkId Output link simplex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeLinkNumber()
      inline int getEdgeLink(const int &edgeId, 
        const int &localLinkId, int &linkId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedEdgeLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeLink query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeLink(
          edgeId, localLinkId, linkId);
      }
      
      /// Get the number of simplicies in the link of the \p edgeId-th edge.
      ///
      /// In 2D, this will return the number of vertices in the link, in 3D the
      /// number of edges.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \return Returns the number of cells in the link of the edge. 
      inline int getEdgeLinkNumber(const int &edgeId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedEdgeLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeLinkNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeLinkNumber(edgeId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link simplices for all edges.
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers representing vertices in 2D and
      /// edges in 3D.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge link list.
      inline const vector<vector<int> > *getEdgeLinks(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedEdgeLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeLinks query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getEdgeLinks();
      }
      
      /// Get the \p localStarId-th cell of the star of the \p edgeId-th edge.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param edgeId Input global edge identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getEdgeStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeStarNumber()
      inline int getEdgeStar(const int &edgeId,
        const int &localStarId, int &starId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedEdgeStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeStar query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeStar(
          edgeId, localStarId, starId);
      }
      
      /// Get the number of star cells for the \p edgeId-th edge.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier
      /// \return Returns the number of star cells.
      inline int getEdgeStarNumber(const int &edgeId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedEdgeStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeStarNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeStarNumber(edgeId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all edges.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding edge.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge star list.
      inline const vector<vector<int> > *getEdgeStars(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedEdgeStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeStars query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getEdgeStars();
      }
      
      /// Get the \p localTriangleId-th triangle id of the \p edgeId-th edge.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param edgeId Input global edge identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getEdgeTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeTriangleNumber()
      inline int getEdgeTriangle(const int &edgeId,
        const int &localTriangleId, int &triangleId) const{
         
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedEdgeTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeTriangle query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeTriangle(
          edgeId, localTriangleId, triangleId);
      }
      
      /// Get the number of triangles for the \p edgeId-th edge.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \return Returns the number of edge triangles.
      inline int getEdgeTriangleNumber(const int &edgeId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
      
        if(!abstractTriangulation_->hasPreprocessedEdgeTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeTriangleNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeTriangleNumber(edgeId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all edges.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding edge.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge triangle list.
      inline const vector<vector<int> > *getEdgeTriangles(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedEdgeTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeTriangles query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getEdgeTriangles();
      }
      
      /// Get the \p localVertexId-th vertex identifier of the \p edgeId-th 
      /// edge.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \param localVertexId Input local vertex identifier (0 or 1).
      /// \param vertexId Output global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getEdgeVertex(const int &edgeId, 
        const int &localVertexId, int &vertexId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "EdgeVertex query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getEdgeVertex(
          edgeId, localVertexId, vertexId);
      }
      
      /// Get the dimensions of the grid if the current object is the implicit
      /// triangulation of a regular grid.
      /// \param dimensions Vector that will be filled with the dimensions of 
      /// the grid. This vector has 3 entries (first: x, second: y, third: z).
      /// \return Returns 0 upon success, negative values otherwise (for 
      /// instance, if the object is not representing a regular grid).
      inline int getGridDimensions(vector<int> &dimensions){
        
        if((gridDimensions_[0] == -1)
          &&(gridDimensions_[1] == -1)
          &&(gridDimensions_[2] == -1)){
          return -1;
        }
        
        dimensions.resize(3);
        dimensions[0] = gridDimensions_[0];
        dimensions[1] = gridDimensions_[1];
        dimensions[2] = gridDimensions_[2];
        
        return 0;
      }
      
      /// Get the number of cells in the triangulation.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \return Returns the number of cells.
      inline int getNumberOfCells() const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->getNumberOfCells();
      }
      
      /// Get the number of edges in the triangulation.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns the number of edges.
      inline int getNumberOfEdges() const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "NumberOfEdges query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getNumberOfEdges();
      }
      
      /// Get the number of triangles in the triangulation.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns the number of triangles.
      inline int getNumberOfTriangles() const{
#ifndef withKamikaze
        
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "NumberOfTriangles query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getNumberOfTriangles();
      }
      
      /// Get the number of vertices in the triangulation.
      /// \return Returns the number of vertices.
      inline int getNumberOfVertices() const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        return abstractTriangulation_->getNumberOfVertices();
      }
     
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles of the triangulation.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of vertex identifiers.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle list.
      inline const vector<vector<int> > *getTriangles(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "Triangles query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getTriangles();
      }
      
      /// Get the \p localEdgeId-th edge of the \p triangleId-th triangle.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \param localEdgeId Input local edge identifier, 
      /// in [0, getTriangleEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleEdgeNumber()
      inline int getTriangleEdge(const int &triangleId, 
        const int &localEdgeId, int &edgeId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedTriangleEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleEdge query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleEdge(
          triangleId, localEdgeId, edgeId);
      }
      
      /// Get the number of edges of the \p triangleId-th triangle.
      ///
      /// Here, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of cells in the link of the triangle. 
      inline int getTriangleEdgeNumber(const int &triangleId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleEdgeNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleEdgeNumber(triangleId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges for all triangles.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles. Each entry is a vector of identifiers representing the
      /// edges connected to the triangle (3).
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle edge list.
      inline const vector<vector<int> > *getTriangleEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleEdges query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getTriangleEdges();
      }
      
      /// Get the \p localLinkId-th simplex of the link of the \p triangleId-th 
      /// triangle.
      ///
      /// Here, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// In 3D, the output \p linkId refers to a vertex identifier.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getTriangleLinkNumber()].
      /// \param linkId Output link simplex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleLinkNumber()
      inline int getTriangleLink(const int &triangleId, 
        const int &localLinkId, int &linkId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedTriangleLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleLink query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -1;
        }
#endif
        return abstractTriangulation_->getTriangleLink(
          triangleId, localLinkId, linkId);
      }
      
      /// Get the number of simplices in the link of the \p triangleId-th 
      /// triangle.
      ///
      /// Here, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// In 3D, this will return the number of vertices in the link.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of simplices in the link of the triangle. 
      inline int getTriangleLinkNumber(const int &triangleId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleLinkNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleLinkNumber(triangleId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link simplices for all triangles.
      ///
      /// Here, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of identifiers representing a vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle link list.
      inline const vector<vector<int> > *getTriangleLinks(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleLinks query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getTriangleLinks();
      }
     
      /// Get the \p localStarId-th cell of the star of the \p triangleId-th 
      /// triangle.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation has
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param triangleId Input global triangle identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getTriangleStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleStarNumber()
      inline int getTriangleStar(const int &triangleId,
        const int &localStarId, int &starId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedTriangleStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleStar query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleStar(
          triangleId, localStarId, starId);
      }
     
      /// Get the number of star cells for the \p triangleId-th triangle.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also,the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of star cells.
      inline int getTriangleStarNumber(const int &triangleId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleStarNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleStarNumber(triangleId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all triangles.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding triangle.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle star list.
      inline const vector<vector<int> > *getTriangleStars(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedTriangleStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleStars query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getTriangleStars();
      }
     
      /// Get the \p localVertexId-th vertex identifier of the \p triangleId-th 
      /// triangle.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global edge identifier.
      /// \param localVertexId Input local vertex identifier (in [0, 2]).
      /// \param vertexId Output global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getTriangleVertex(const int &triangleId,
        const int &localVertexId, int &vertexId) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "TriangleVertex query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getTriangleVertex(
          triangleId, localVertexId, vertexId);
      }
     
      /// Get the \p localEdgeId-th edge identifier connected to the 
      /// \p vertexId-th 
      /// vertex.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localEdgeId Input local edge identifier,
      /// in [0, getVertexEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexEdgeNumber()
      inline int getVertexEdge(const int &vertexId, 
        const int &localEdgeId, int &edgeId) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedVertexEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexEdge query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexEdge(
          vertexId, localEdgeId, edgeId);
      }
     
      /// Get the number of edges connected to the \p vertexId-th vertex.
      ///
      /// Here,the notion of edge only makes sense if the triangulation has 
      /// a dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of edges connected to the vertex.
      inline int getVertexEdgeNumber(const int &vertexId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedVertexEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexEdgeNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexEdgeNumber(vertexId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edge identifiers for all vertices.
      ///
      /// Here, the notion of edge only makes sense if the triangulation
      /// has a dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of edges connected to the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex edge list.
      inline const vector<vector<int> > *getVertexEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedVertexEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexEdges query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getVertexEdges();
      }
      
      /// Get the \p localLinkId-th simplex of the link of the \p vertexId-th 
      /// vertex.
      ///
      /// The output \p linkId refers in 2D to an edge identifier and in 3D to 
      /// a triangle identifier.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getVertexLinkNumber()].
      /// \param linkId Output link simplex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexLinkNumber()
      inline int getVertexLink(const int &vertexId, 
        const int &localLinkId, int &linkId) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedVertexLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexLink query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexLink(
          vertexId, localLinkId, linkId);
      }
      
      /// Get the number of simplices in the link of the \p vertexId-th vertex.
      ///
      /// In 2D, this will return the number of edges in the link, in 3D the 
      /// number of triangles.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of cells in the link of the vertex. 
      inline int getVertexLinkNumber(const int &vertexId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedVertexLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexLinkNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexLinkNumber(vertexId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link simplices for all vertices.
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers representing edges in 2D and 
      /// triangles in 3D.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex link list.
      inline const vector<vector<int> > *getVertexLinks(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedVertexLinks()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexLinks query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getVertexLinks();
      }
      
      /// Get the \p localNeighborId-th vertex neighbor of the \p vertexId-th 
      /// vertex.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localNeighborId Input local neighbor identifier,
      /// in [0, getVertexNeighborNumber()].
      /// \param neighborId Output global neighbor vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexNeighborNumber()
      inline int getVertexNeighbor(const int &vertexId, 
        const int &localNeighborId, int &neighborId) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexNeighbor query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexNeighbor(
          vertexId, localNeighborId, neighborId);
      }
      
      /// Get the number of vertex neighbors for the \p vertexId-th vertex.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number vertex neighbors.
      inline int getVertexNeighborNumber(const int &vertexId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexNeighborNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif     
        return abstractTriangulation_->getVertexNeighborNumber(vertexId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of vertex neighbor identifiers for all vertices.
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of vertex neighbors for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex neighbor list.
      inline const vector<vector<int> > *getVertexNeighbors(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexNeighbors query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif 
        return abstractTriangulation_->getVertexNeighbors();
      }
      
      /// Get the point (3D coordinates) for the \p vertexId-th vertex.
      /// \param vertexId Input global vertex identifier.
      /// \param x Output x coordinate.
      /// \param y Output y coordinate.
      /// \param z Output z coordinate.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getVertexPoint(const int &vertexId, 
        float &x, float &y, float &z) const{
          
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        return abstractTriangulation_->getVertexPoint(vertexId, x, y, z);
      }
      
      /// Get the \p localStarId-th cell of the star of the \p vertexId-th 
      /// vertex.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param vertexId Input global vertex identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getVertexStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexStarNumber()
      inline int getVertexStar(const int &vertexId, 
        const int &localStarId, int &starId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedVertexStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexStar query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexStar() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif          
        return abstractTriangulation_->getVertexStar(
          vertexId, localStarId, starId);
      }
        
      /// Get the number of star cells for the \p vertexId-th vertex.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier
      /// \return Returns the number of star cells.
      inline int getVertexStarNumber(const int &vertexId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
      
        if(!abstractTriangulation_->hasPreprocessedVertexStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexStarNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexStarNumber(vertexId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all vertices.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex star list.
      inline const vector<vector<int> > *getVertexStars(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedVertexStars()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexStars query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexStars() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif         
        return abstractTriangulation_->getVertexStars();
      }
      
      /// Get the \p localTriangleId-th triangle id of the 
      /// \p vertexId-th vertex.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param vertexId Input global vertex identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getVertexTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexTriangleNumber()
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getVertexTriangle(const int &vertexId,
        const int &localTriangleId, int &triangleId) const{
         
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
          
        if(!abstractTriangulation_->hasPreprocessedVertexTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexTriangle query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexTriangle(
          vertexId, localTriangleId, triangleId);
      }
      
      /// Get the number of triangles for the \p vertexId-th vertex.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of vertex triangles.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getVertexTriangleNumber(const int &vertexId) const{
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
        
        if(!abstractTriangulation_->hasPreprocessedVertexTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexTriangleNumber query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return -2;
        }
#endif
        return abstractTriangulation_->getVertexTriangleNumber(vertexId);
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all vertices.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex triangle list.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline const vector<vector<int> > *getVertexTriangles(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return NULL;
        
        if(!abstractTriangulation_->hasPreprocessedVertexTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "VertexTriangles query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return NULL;
        }
#endif
        return abstractTriangulation_->getVertexTriangles();
      }
      
      /// Check if the edge with global identifier \p edgeId is on the boundary 
      /// of the domain. 
      /// 
      /// For 2D triangulations, this function will return true if the edge is 
      /// a boundary edge. For 3D triangulations, this function will return 
      /// true if the edge belongs to a boundary triangle.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessBoundaryEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \return Returns true if the edge is on the boundary, false otherwise.
      inline bool isEdgeOnBoundary(const int &edgeId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return false;
        
        if(!abstractTriangulation_->hasPreprocessedBoundaryEdges()){
          stringstream msg;
          msg << "[Triangulation] "
            << "BoundaryEdge query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessBoundaryEdges() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return false;
        }
#endif
        return abstractTriangulation_->isEdgeOnBoundary(edgeId);
      }
      
      /// Check if the data structure is empty or not.
      /// \return Returns true if empty, false otherwise.
      inline bool isEmpty() const {
        return !abstractTriangulation_;
      }
      
      /// Check if the triangle with global identifier \p triangleId is on the 
      /// boundary of the domain. 
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// For 2D triangulations, this function will return false all the time.
      /// For 3D triangulations, this function will return true if the triangle
      /// is a boundary triangle.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessBoundaryTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns true if the triangle is on the boundary, false 
      /// otherwise.
      inline bool isTriangleOnBoundary(const int &triangleId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return false;
        
        if(!abstractTriangulation_->hasPreprocessedBoundaryTriangles()){
          stringstream msg;
          msg << "[Triangulation] "
            << "BoundaryTriangle query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessBoundaryTriangles() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return false;
        }
#endif
        return abstractTriangulation_->isTriangleOnBoundary(triangleId);
      }
      
      /// Check if the vertex with global identifier \p vertexId is on the 
      /// boundary of the domain. 
      /// 
      /// For 2D triangulations, this function will return true if the vertex 
      /// belongs to a boundary edge. For 3D triangulations, this function will 
      /// return true if the vertex belongs to a boundary triangle.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessBoundaryVertices() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns true if the vertex is on the boundary, false 
      /// otherwise.
      inline bool isVertexOnBoundary(const int &vertexId) const{
#ifndef withKamikaze
        if(isEmptyCheck())
          return false;
        
        if(!abstractTriangulation_->hasPreprocessedBoundaryVertices()){
          stringstream msg;
          msg << "[Triangulation] "
            << "BoundaryVertex query without pre-process!"
            << endl;
          msg << "[Triangulation] "
            << "Please call preprocessBoundaryVertices() in a"
            << " pre-process." << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          return false;
        }
#endif
        return abstractTriangulation_->isVertexOnBoundary(vertexId);
      }
      
      /// Pre-process the boundary edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following function(s):
      ///   - isEdgeOnBoundary()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa isEdgeOnBoundary()
      inline int preprocessBoundaryEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessBoundaryEdges();
      }
      
      /// Pre-process the boundary triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following function(s):
      ///   - isTriangleOnBoundary()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa isTriangleOnBoundary()
      inline int preprocessBoundaryTriangles(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessBoundaryTriangles();
      }
      
      /// Pre-process the boundary vertices.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following function(s):
      ///   - isVertexOnBoundary()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa isVertexOnBoundary()
      inline int preprocessBoundaryVertices(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessBoundaryVertices();
      }
     
      /// Pre-process the cell edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellEdge()
      ///   - getCellEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellEdge()
      /// \sa getCellEdgeNumber()
      inline int preprocessCellEdges(){
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessCellEdges();
      }

      /// Pre-process the cell neighbors.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellNeighbor()
      ///   - getCellNeighbors()
      ///   - getCellNeighborNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellNeighbor()
      /// \sa getCellNeighbors()
      /// \sa getCellNeighborNumber()
      inline int preprocessCellNeighbors(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessCellNeighbors();
      }
      
      /// Pre-process the cell triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellTriangle()
      ///   - getCellTriangles()
      ///   - getCellTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellTriangle()
      /// \sa getCellTriangles()
      /// \sa getCellTriangleNumber()
      inline int preprocessCellTriangles(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif        
        
        return abstractTriangulation_->preprocessCellTriangles();
      }
      
      /// Pre-process the edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdges()
      ///   - getEdgeVertex()
      ///   - getNumberOfEdges()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdges()
      /// \sa getEdgeVertex()
      /// \sa getNumberOfEdges()
      inline int preprocessEdges(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif
        
        return abstractTriangulation_->preprocessEdges();
      }
      
      /// Pre-process the edge links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeLink()
      ///   - getEdgeLinks()
      ///   - getEdgeLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeLink()
      /// \sa getEdgeLinks()
      /// \sa getEdgeLinkNumber()
      inline int preprocessEdgeLinks(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif        
        
        return abstractTriangulation_->preprocessEdgeLinks();
      }
      
      /// Pre-process the edge stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeStar()
      ///   - getEdgeStars()
      ///   - getEdgeStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeStar()
      /// \sa getEdgeStars()
      /// \sa getEdgeStarNumber()
      inline int preprocessEdgeStars(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif          
        
        return abstractTriangulation_->preprocessEdgeStars();
      }
      
      /// Pre-process the edge triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeTriangle()
      ///   - getEdgeTriangles()
      ///   - getEdgeTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeTriangle()
      /// \sa getEdgeTriangles()
      /// \sa getEdgeTriangleNumber()
      inline int preprocessEdgeTriangles(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif         
        
        return abstractTriangulation_->preprocessEdgeTriangles();
      }
      
      /// Pre-process the triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getNumberOfTriangles()
      ///   - getTriangles()
      ///   - getTriangleVertex()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getNumberOfTriangles()
      /// \sa getTriangles()
      /// \sa getTriangleVertex()
      inline int preprocessTriangles(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif    
        
        return abstractTriangulation_->preprocessTriangles();
      }
      
      /// Pre-process the triangle edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleEdge()
      ///   - getTriangleEdges()
      ///   - getTriangleEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleEdge()
      /// \sa getTriangleEdges()
      /// \sa getTriangleEdgeNumber()
      inline int preprocessTriangleEdges(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif 
        
        return abstractTriangulation_->preprocessTriangleEdges();
      }
      
      /// Pre-process the triangle links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleLink()
      ///   - getTriangleLinks()
      ///   - getTriangleLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleLink()
      /// \sa getTriangleLinks()
      /// \sa getTriangleLinkNumber()
      inline int preprocessTriangleLinks(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif 
        
        return abstractTriangulation_->preprocessTriangleLinks();
      }
      
      /// Pre-process the triangle stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleStar()
      ///   - getTriangleStars()
      ///   - getTriangleStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleStar()
      /// \sa getTriangleStars()
      /// \sa getTriangleStarNumber()
      inline int preprocessTriangleStars(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif         
        
        return abstractTriangulation_->preprocessTriangleStars();
      }
      
      /// Pre-process the vertex edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexEdge()
      ///   - getVertexEdges()
      ///   - getVertexEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexEdge()
      /// \sa getVertexEdges()
      /// \sa getVertexEdgeNumber()
      inline int preprocessVertexEdges(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif         
        
        return abstractTriangulation_->preprocessVertexEdges();
      }
      
      /// Pre-process the vertex links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexLink()
      ///   - getVertexLinks()
      ///   - getVertexLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexLink()
      /// \sa getVertexLinks()
      /// \sa getVertexLinkNumber()
      inline int preprocessVertexLinks(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif         
        
        return abstractTriangulation_->preprocessVertexLinks();
      }
      
      /// Pre-process the vertex neighbors.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexNeighbor()
      ///   - getVertexNeighbors()
      ///   - getVertexNeighborNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexNeighbor()
      /// \sa getVertexNeighbors()
      /// \sa getVertexNeighborNumber()
      inline int preprocessVertexNeighbors(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif         
        
        return abstractTriangulation_->preprocessVertexNeighbors();
      }
      
      /// Pre-process the vertex stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexStar()
      ///   - getVertexStars()
      ///   - getVertexStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexStar()
      /// \sa getVertexStars()
      /// \sa getVertexStarNumber()
      inline int preprocessVertexStars(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif          
        
        return abstractTriangulation_->preprocessVertexStars();
      }
      
      /// Pre-process the vertex triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexTriangle()
      ///   - getVertexTriangles()
      ///   - getVertexTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexTriangle()
      /// \sa getVertexTriangles()
      /// \sa getVertexTriangleNumber()
      inline int preprocessVertexTriangles(){
        
#ifndef withKamikaze
        if(isEmptyCheck())
          return -1;
#endif        
        
        return abstractTriangulation_->preprocessVertexTriangles();
      }
     
      /// Tune the debug level (default: 0)
      inline int setDebugLevel(const int &debugLevel){
        explicitTriangulation_.setDebugLevel(debugLevel);
        implicitTriangulation_.setDebugLevel(debugLevel);
        debugLevel_ = debugLevel;
        return 0;
      }

      /// Set the input cells for the triangulation.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \param cellNumber Number of input cells.
      /// \param cellArray Pointer to the input cells. This pointer should point
      /// to an array of long long int where cells are stored one after the 
      /// other. In particular, each cell starts by the number of vertices in 
      /// it, followed by the identifiers of its vertices. This corresponds to 
      /// the default cell array representation in VTK.
      /// \return Returns 0 upon success, negative values otherwise.
      ///
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      ///
      /// \warning If this ttk::Triangulation object is already representing a 
      /// valid triangulation, this information will be over-written (which 
      /// means that pre-processing functions should be called again).
      inline int setInputCells(const int &cellNumber,
        const long long int *cellArray){
        
        abstractTriangulation_ = &explicitTriangulation_;
        gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
        
        return explicitTriangulation_.setInputCells(cellNumber, cellArray);
      }
      
      /// Set the specifications of the input grid to implicitly represent as a
      /// triangulation.
      /// \param xOrigin Input x coordinate of the grid origin.
      /// \param yOrigin Input y coordinate of the grid origin.
      /// \param zOrigin Input z coordinate of the grid origin.
      /// \param xSpacing Input spacing along the x dimension.
      /// \param ySpacing Input spacing along the y dimension.
      /// \param zSpacing Input spacing along the z dimension.
      /// \param xDim Input number of vertices along the x dimension.
      /// \param yDim Input number of vertices along the y dimension.
      /// \param zDim Input number of vertices along the z dimension.
      /// \return Returns 0 upon success, negative values otherwise.
      ///
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      ///
      /// \warning If this ttk::Triangulation object is already representing a 
      /// valid triangulation, this information will be over-written (which 
      /// means that pre-processing functions should be called again).
      inline int setInputGrid(
        const float &xOrigin, const float &yOrigin, const float &zOrigin,
        const float &xSpacing, const float &ySpacing, const float &zSpacing,
        const int &xDim, const int &yDim, const int &zDim){
      
        abstractTriangulation_ = &implicitTriangulation_;
        
        gridDimensions_[0] = xDim;
        gridDimensions_[1] = yDim;
        gridDimensions_[2] = zDim;
        
        return implicitTriangulation_.setInputGrid(
          xOrigin, yOrigin, zOrigin,
          xSpacing, ySpacing, zSpacing,
          xDim, yDim, zDim);
        return 0;
      }
      
      /// Set the input 3D points of the triangulation.
      /// \param pointNumber Number of input vertices.
      /// \param pointSet Pointer to the 3D points. This pointer should point to
      /// an array of float where points are stored one after the other. 
      /// In particular, each point is represented by X-Y-Z coordinates (one 
      /// after the other). This corresponds to the default point set 
      /// representation in VTK.
      /// \return Returns 0 upon success, negative values otherwise.
      ///
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      ///
      /// \warning If this ttk::Triangulation object is already representing a 
      /// valid triangulation, this information will be over-written (which 
      /// means that pre-processing functions should be called again).
      inline int setInputPoints(const int &pointNumber, const float *pointSet){
        
        abstractTriangulation_ = &explicitTriangulation_;
        gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
        return explicitTriangulation_.setInputPoints(pointNumber, pointSet);
      }

      /// Tune the number of active threads (default: number of logical cores)
      inline int setThreadNumber(const int &threadNumber){
        explicitTriangulation_.setThreadNumber(threadNumber);
        implicitTriangulation_.setThreadNumber(threadNumber);
        threadNumber_ = threadNumber;
        return 0;
      }

      /// Internal usage. Pass the execution context (debug level, number of 
      /// threads, etc.) to the implementing classes.
      inline int setWrapper(const Wrapper *wrapper){
        explicitTriangulation_.setWrapper(wrapper);
        implicitTriangulation_.setWrapper(wrapper);
        return 0;
      }
      
    protected:
   
      inline bool isEmptyCheck() const{
        if(!abstractTriangulation_){
          stringstream msg;
          msg << "[Triangulation] Trying to access an empty data-structure!"
            << endl;
          dMsg(cerr, msg.str(), fatalMsg);
          return true;
        }
        return false;
      }
      
      int                 gridDimensions_[3];
      
      AbstractTriangulation
                          *abstractTriangulation_;
      ExplicitTriangulation 
                          explicitTriangulation_;
      ImplicitTriangulation
                          implicitTriangulation_;
  };
}

// if the package is not a template, comment the following line
// #include                  <Triangulation.cpp>

#endif // _TRIANGULATION_H
