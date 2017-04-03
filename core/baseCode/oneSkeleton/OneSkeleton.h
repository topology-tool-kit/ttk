/// \ingroup baseCode
/// \class ttk::OneSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief %OneSkeleton processing package.
///
/// %OneSkeleton is a processing package that handles the 1-skeleton (edges)
/// of a triangulation.
/// \sa Triangulation
/// \sa vtkTriangulation
/// \sa vtkOneSkeleton

#ifndef _ONESKELETON_H
#define _ONESKELETON_H

#include                  <map>

// base code includes
#include                  <ZeroSkeleton.h>
#include                  <Wrapper.h>

namespace ttk{
  
  class OneSkeleton : public Debug{

    public:
        
      OneSkeleton();
      
      ~OneSkeleton();

      /// Compute the link of each edge of a 2D triangulation (unspecified 
      /// behavior if the input mesh is not a valid triangulation).
      /// \param edgeList List of edges. The size of this vector 
      /// should be equal to the number of edges in the triangulation. Each
      /// entry is a pair of vertex identifiers.
      /// \param edgeStars List of edge stars. The size of this vector should be
      /// equal to the number of edges. Each entry is a vector of triangle
      /// identifiers.
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param edgeLinks Output edge links. The size of this vector 
      /// will be equal to the number of edges in the triangulation. Each 
      /// entry will be a vector listing the vertices in the link of the 
      /// corresponding vertex.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeLinks(const vector<pair<int, int> > &edgeList,
        const vector<vector<int> > &edgeStars,
        const long long int *cellArray,
        vector<vector<int> > &edgeLinks) const;
      
      /// Compute the link of each edge of a 3D triangulation (unspecified 
      /// behavior if the input mesh is not a valid triangulation).
      /// \param edgeList List of edges. The size of this vector 
      /// should be equal to the number of edges in the triangulation. Each
      /// entry is a pair of vertex identifiers.
      /// \param edgeStars List of edge stars. The size of this vector should be
      /// equal to the number of edges. Each entry is a vector of tetrahedron
      /// identifiers.
      /// \param cellEdges List of celle edges. The size of this vector 
      /// should be equal to the number of tetrahedra in the triangulation. Each
      /// entry is a vector of edge identifiers.
      /// \param edgeLinks Output edge links. The size of this vector 
      /// will be equal to the number of edges in the triangulation. Each 
      /// entry will be a vector listing the vertices in the link of the 
      /// corresponding vertex.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeLinks(const vector<pair<int, int> > &edgeList,
        const vector<vector<int> > &edgeStars,
        const vector<vector<int> > &cellEdges,
        vector<vector<int> > &edgeLinks) const;
      
      /// Compute the list of edges of a valid triangulation.
      /// \param vertexNumber Number of vertices in the triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param edgeList Output edge list (each entry is an ordered pair of 
      /// vertex identifiers).
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeList(const int &vertexNumber, const int &cellNumber, 
        const long long int *cellArray,
        vector<pair<int, int> > &edgeList) const;
      
      /// Compute the list of edges of multiple triangulations.
      /// \param cellArrays Vector of cells. For each triangulation, each entry
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param edgeList Output edge list (each entry is an ordered pair of 
      /// vertex identifiers).
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeLists(const vector<vector<long long int> >  &cellArrays,
        vector<vector<pair<int, int> > > &edgeLists) const;
      
          
      /// Compute the 3-star of all the edges of a triangulation (for each
      /// edge, list of the 3-dimensional cells connected to it).
      /// \param vertexNumber Number of vertices in the triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param starList Output list of 3-stars. The size of this vector will
      /// be equal to the number of edges in the mesh. Each entry stores a 
      /// vector that lists the identifiers of all 3-dimensional cells 
      /// connected to the entry's edge.
      /// \param edgeList Optional list of edges. If NULL, the function will 
      /// compute this list anyway and free the related memory upon return.
      /// If not NULL but pointing to an empty vector, the function will fill 
      /// this empty vector (useful if this list needs to be used later on by 
      /// the calling program). If not NULL but pointing to a non-empty vector,
      /// this function will use this vector as internal edge list. If this 
      /// vector is not empty but incorrect, the behavior is unspecified.
      /// \param vertexStars Optional list of vertex stars (list of 
      /// 3-dimensional cells connected to each vertex). If NULL, the 
      /// function will compute this list anyway and free the related memory
      /// upon return. If not NULL but pointing to an empty vector, the 
      /// function will fill this empty vector (useful if this list needs 
      /// to be used later on by the calling program). If not NULL but pointing
      /// to a non-empty vector, this function will use this vector as internal 
      /// vertex star list. If this vector is not empty but incorrect, the 
      /// behavior is unspecified.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeStars(const int &vertexNumber, const int &cellNumber,
        const long long int *cellArray,
        vector<vector<int> > &starList,
        vector<pair<int, int> > *edgeList = NULL,
        vector<vector<int> > *vertexStars = NULL) const;
      
      /// Compute the list of edges of a sub-portion of a valid triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// considered subset of the triangulation (number of tetrahedra in 3D, 
      /// triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param edgeList Output edge list (each entry is an ordered pair of 
      /// vertex identifiers).
      /// \return Returns 0 upon success, negative values otherwise.
      int buildEdgeSubList(
        const int &cellNumber, const long long int *cellArray,
        vector<pair<int, int> > &edgeList) const;
        
      
    protected:
    
  };
}

// if the package is not a template, comment the following line
// #include                  <OneSkeleton.cpp>

#endif // ONESKELETON_H
