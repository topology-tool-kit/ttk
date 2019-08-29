/// \ingroup base
/// \class ttk::ManifoldCheck
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2018.
///
/// \brief TTK processing package for manifold checks.
///
/// This class performs a manifold check for each simplex, by counting the
/// number of connected components of link. On a d-dimensional triangulation,
/// this number should be equal to 1 for all but (d-1)-simplices, for which it
/// can be 1 (boundary simplices) or 2 (interior simplices). Other values
/// indicate a non-manifold simplex.
///
/// The link component number is stored as an integer array for each type of
/// simplex.
///
/// \sa ttk::Triangulation
/// \sa ttkManifoldCheck.cpp %for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <UnionFind.h>
#include <Wrapper.h>

namespace ttk {

  class ManifoldCheck : public Debug {

  public:
    ManifoldCheck();

    ~ManifoldCheck();

    /// Execute the package.
    int execute() const;

    /// Register the output vector for vertex link component number
    inline int setVertexLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *vertexVector) {
      vertexLinkComponentNumber_ = vertexVector;
      return 0;
    }

    /// Register the output vector for edge link component number
    inline int setEdgeLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *edgeVector) {
      edgeLinkComponentNumber_ = edgeVector;
      return 0;
    }

    /// Register the output vector for triangle link component number
    inline int setTriangleLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *triangleVector) {
      triangleLinkComponentNumber_ = triangleVector;
      return 0;
    }

    /// Setup a (valid) triangulation object for this TTK base object.
    ///
    /// \pre This function should be called prior to any usage of this TTK
    /// object, in a clearly distinct pre-processing step that involves no
    /// traversal or computation at all. An error will be returned otherwise.
    ///
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement. Therefore, it is recommended to
    /// call this function ONLY in the pre-processing steps of your program.
    /// Note however, that your triangulation object must be valid when
    /// calling this function (i.e. you should have filled it at this point,
    /// see the setInput*() functions of ttk::Triangulation).
    /// See ttkManifoldCheck
    /// for further examples.
    ///
    /// \param triangulation Pointer to a valid triangulation.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa ttk::Triangulation
    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;

      if(triangulation_) {

        triangulation_->preprocessVertexLinks();
        triangulation_->preprocessEdgeLinks();
        triangulation_->preprocessTriangleLinks();
      }

      return 0;
    }

  protected:
    int vertexManifoldCheck(const ttk::SimplexId &vertexId) const;

    int edgeManifoldCheck(const ttk::SimplexId &edgeId) const;

    std::vector<ttk::SimplexId> *vertexLinkComponentNumber_;
    std::vector<ttk::SimplexId> *edgeLinkComponentNumber_;
    std::vector<ttk::SimplexId> *triangleLinkComponentNumber_;
    Triangulation *triangulation_;
  };
} // namespace ttk
