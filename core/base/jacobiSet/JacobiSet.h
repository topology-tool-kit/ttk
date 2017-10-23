/// \ingroup baseCode
/// \class ttk::JacobiSet 
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK processing package for the computation of the Jacobi set of 
/// bivariate volumetric data.
///
/// Given a bivariate scalar field defined on a PL 3-manifold, this package
/// produces the list of Jacobi edges (each entry is a pair given by the edge
/// identifier and the Jacobi edge type).
/// \param dataTypeU Data type of the input first component field (char, float, 
/// etc.).
/// \param dataTypeV Data type of the input second component field (char, float,
/// etc.)
/// 
/// \b Related \b publication \n
/// "Jacobi sets of multiple Morse functions" \n
/// Herbert Edelsbrunner, John Harer \n
/// Foundations of Computational Mathematics. Cambridge University Press, 2002.
///
/// \sa vtkJacobiSet.cpp %for a usage example.

#ifndef _JACOBISET_H
#define _JACOBISET_H

// base code includes
#include                  <ScalarFieldCriticalPoints.h>
#include                  <Triangulation.h>
#include                  <UnionFind.h>
#include                  <Wrapper.h>


namespace ttk{
  
  template <class dataTypeU, class dataTypeV> class JacobiSet : public Debug{

    public:
        
      JacobiSet();
      
      ~JacobiSet();

      int connectivityPreprocessing(const vector<vector<int> > &edgeStarList,
        vector<vector<pair<int, int> > > &edgeFanLinkEdgeLists,
        vector<vector<long long int> > &edgeFans,
        vector<int> &sosOffsets) const;
      
      int execute(vector<pair<int, char> > &jacobiSet);
    
      char getCriticalType(const int &edgeId);
      
      int perturbate(const dataTypeU &uEpsilon = pow(10, -DBL_DIG),
        const dataTypeV &vEpsilon = pow(10, -DBL_DIG)) const;
      
      int setEdgeFans(const vector<vector<long long int> > *edgeFans){
        edgeFans_ = edgeFans;
        return 0;
      }
      
      int setEdgeFanLinkEdgeList(
        const vector<vector<pair<int, int> > > *edgeFanLinkEdgeLists){
        edgeFanLinkEdgeLists_ = edgeFanLinkEdgeLists;
        return 0;
      }
      
      int setEdgeList(const vector<pair<int, int> > *edgeList){
        edgeList_ = edgeList;
        return 0;
      }
      
      int setInputField(const void *uField, const void *vField){
        
        uField_ = uField;
        vField_ = vField;
        return 0;
      }
      
      int setSosOffsets(vector<int> *sosOffsets){
        // legacy API
        return setSosOffsetsU(sosOffsets);
      }
      
      int setSosOffsetsU(vector<int> *sosOffsets){
        sosOffsetsU_ = sosOffsets;
        return 0;
      }
      
      int setSosOffsetsV(vector<int> *sosOffsets){
        sosOffsetsV_ = sosOffsets;
        return 0;
      }

      // NOTE: here it's not clear how vtk builds vtkIdType 
      // to check on bigger data-sets
      int setTetList(const long long int *tetList){
        tetList_ = tetList;
        return 0;
      }
      
      int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }
      
      int setupTriangulation(Triangulation *triangulation){
        
        triangulation_ = triangulation;
        
        // pre-condition functions
        if(triangulation_){
          triangulation_->preprocessEdges();
          triangulation_->preprocessEdgeStars();
        }
        
        return 0;
      }
      
    protected:
    
      int executeLegacy(vector<pair<int, char> > &jacobiSet);
      
      int                   vertexNumber_;
      const long long int   *tetList_;
      const void            *uField_, *vField_;
      const vector<pair<int, int> > *edgeList_;
      // for each edge, one skeleton of its triangle fan
      const vector<vector<pair<int, int> > > *edgeFanLinkEdgeLists_;
      // for each edge, the one skeleton of its triangle fan
      const vector<vector<long long int> > *edgeFans_;
      vector<int>           *sosOffsetsU_, *sosOffsetsV_;
      vector<int>           localSosOffsetsU_, localSosOffsetsV_;
      Triangulation         *triangulation_;
  };
}

// if the package is not a template, comment the following line
#include                  <JacobiSet.inl>

#endif // JACOBISET_H
