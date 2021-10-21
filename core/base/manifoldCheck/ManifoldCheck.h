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

namespace ttk {

  class ManifoldCheck : virtual public Debug {

  public:
    ManifoldCheck();

    ~ManifoldCheck();

    /// Execute the package.
    template <class triangulationType = AbstractTriangulation>
    int execute(const triangulationType *triangulation) const;

    /// Register the output std::vector for vertex link component number
    inline int setVertexLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *vertexVector) {
      vertexLinkComponentNumber_ = vertexVector;
      return 0;
    }

    /// Register the output std::vector for edge link component number
    inline int setEdgeLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *edgeVector) {
      edgeLinkComponentNumber_ = edgeVector;
      return 0;
    }

    /// Register the output std::vector for triangle link component number
    inline int setTriangleLinkComponentNumberVector(
      std::vector<ttk::SimplexId> *triangleVector) {
      triangleLinkComponentNumber_ = triangleVector;
      return 0;
    }

    /// Precondition a (valid) triangulation object for this TTK base object.
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
    inline int
      preconditionTriangulation(AbstractTriangulation *const triangulation) {

      if(triangulation) {

        triangulation->preconditionVertexLinks();
        triangulation->preconditionEdgeLinks();
        triangulation->preconditionTriangleLinks();
        triangulation->preconditionVertexEdges();
        triangulation->preconditionVertexTriangles();
      }

      return 0;
    }

  protected:
    template <class triangulationType = AbstractTriangulation>
    int vertexManifoldCheck(const triangulationType *triangulation,
                            const ttk::SimplexId &vertexId) const;

    template <class triangulationType = AbstractTriangulation>
    int edgeManifoldCheck(const triangulationType *triangulation,
                          const ttk::SimplexId &edgeId) const;

    std::vector<ttk::SimplexId> *vertexLinkComponentNumber_;
    std::vector<ttk::SimplexId> *edgeLinkComponentNumber_;
    std::vector<ttk::SimplexId> *triangleLinkComponentNumber_;
  };
} // namespace ttk

template <class triangulationType>
int ttk::ManifoldCheck::execute(const triangulationType *triangulation) const {

  printMsg(ttk::debug::Separator::L1);

  Timer t;

// check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation)
    return -1;
#endif

  SimplexId vertexNumber = triangulation->getNumberOfVertices();

  if(vertexLinkComponentNumber_) {

    vertexLinkComponentNumber_->resize(vertexNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < vertexNumber; i++) {
      (*vertexLinkComponentNumber_)[i] = vertexManifoldCheck(triangulation, i);
    }
  }

  if((edgeLinkComponentNumber_) && (triangulation->getDimensionality() >= 2)) {

    SimplexId edgeNumber = triangulation->getNumberOfEdges();
    edgeLinkComponentNumber_->resize(edgeNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < edgeNumber; i++) {
      (*edgeLinkComponentNumber_)[i] = edgeManifoldCheck(triangulation, i);
    }
  }

  if((triangleLinkComponentNumber_)
     && (triangulation->getDimensionality() == 3)) {

    SimplexId triangleNumber = triangulation->getNumberOfTriangles();
    triangleLinkComponentNumber_->resize(triangleNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < triangleNumber; i++) {
      (*triangleLinkComponentNumber_)[i]
        = triangulation->getTriangleLinkNumber(i);
    }
  }

  printMsg("Processed " + std::to_string(vertexNumber) + " vertices", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 0;
}

template <class triangulationType>
int ttk::ManifoldCheck::vertexManifoldCheck(
  const triangulationType *triangulation, const SimplexId &vertexId) const {

  SimplexId linkSize = triangulation->getVertexLinkNumber(vertexId);

  if(triangulation->getDimensionality() == 1)
    return linkSize;

  std::vector<SimplexId> linkNeighbors;

  for(SimplexId i = 0; i < linkSize; i++) {
    SimplexId linkId = -1;
    triangulation->getVertexLink(vertexId, i, linkId);

    bool isIn = false;
    SimplexId neighborId = -1;

    if(triangulation->getDimensionality() == 2) {
      triangulation->getEdgeVertex(linkId, 0, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation->getEdgeVertex(linkId, 1, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);
    }
    if(triangulation->getDimensionality() == 3) {
      triangulation->getTriangleVertex(linkId, 0, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation->getTriangleVertex(linkId, 1, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation->getTriangleVertex(linkId, 2, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);
    }
  }

  std::vector<UnionFind> seeds(linkNeighbors.size());
  std::vector<UnionFind *> seedList(linkNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)seeds.size(); i++) {
    seedList[i] = &(seeds[i]);
  }

  for(SimplexId i = 0; i < linkSize; i++) {

    SimplexId linkId = -1;
    triangulation->getVertexLink(vertexId, i, linkId);

    SimplexId neighborId0 = -1, neighborId1 = -1, neighborId2 = -1;
    SimplexId uf0 = -1, uf1 = -1, uf2 = -1;

    if(triangulation->getDimensionality() == 2) {
      triangulation->getEdgeVertex(linkId, 0, neighborId0);
      triangulation->getEdgeVertex(linkId, 1, neighborId1);

      // connect the two uf together
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId0) {
          uf0 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId1) {
          uf1 = j;
          break;
        }
      }

      seedList[uf0] = UnionFind::makeUnion(seedList[uf0], seedList[uf1]);
      seedList[uf1] = seedList[uf0];
    }

    if(triangulation->getDimensionality() == 3) {
      triangulation->getTriangleVertex(linkId, 0, neighborId0);
      triangulation->getTriangleVertex(linkId, 1, neighborId1);
      triangulation->getTriangleVertex(linkId, 2, neighborId2);

      // connect the two uf together
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId0) {
          uf0 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId1) {
          uf1 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId2) {
          uf2 = j;
          break;
        }
      }

      seedList[uf0] = UnionFind::makeUnion(seedList[uf0], seedList[uf1]);
      seedList[uf0] = UnionFind::makeUnion(seedList[uf0], seedList[uf2]);
      seedList[uf1] = seedList[uf0];
      seedList[uf2] = seedList[uf0];
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)seedList.size(); i++) {
    seedList[i] = seedList[i]->find();
  }

  std::vector<UnionFind *>::iterator it;
  sort(seedList.begin(), seedList.end());
  it = unique(seedList.begin(), seedList.end());
  seedList.resize(distance(seedList.begin(), it));

  return (SimplexId)seedList.size();
}

template <class triangulationType>
int ttk::ManifoldCheck::edgeManifoldCheck(
  const triangulationType *triangulation, const SimplexId &edgeId) const {

  SimplexId linkSize = triangulation->getEdgeLinkNumber(edgeId);

  if(triangulation->getDimensionality() == 2)
    return linkSize;

  std::vector<SimplexId> linkNeighbors;

  for(SimplexId i = 0; i < linkSize; i++) {
    SimplexId linkId = -1;
    triangulation->getEdgeLink(edgeId, i, linkId);

    bool isIn = false;
    SimplexId neighborId = -1;

    triangulation->getEdgeVertex(linkId, 0, neighborId);
    isIn = false;
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId) {
        isIn = true;
        break;
      }
    }
    if(!isIn)
      linkNeighbors.push_back(neighborId);

    triangulation->getEdgeVertex(linkId, 1, neighborId);
    isIn = false;
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId) {
        isIn = true;
        break;
      }
    }
    if(!isIn)
      linkNeighbors.push_back(neighborId);
  }

  std::vector<UnionFind> seeds(linkNeighbors.size());
  std::vector<UnionFind *> seedList(linkNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)seeds.size(); i++) {
    seedList[i] = &(seeds[i]);
  }

  for(SimplexId i = 0; i < linkSize; i++) {

    SimplexId linkId = -1;
    triangulation->getEdgeLink(edgeId, i, linkId);

    SimplexId neighborId0 = -1, neighborId1 = -1;
    SimplexId uf0 = -1, uf1 = -1;

    triangulation->getEdgeVertex(linkId, 0, neighborId0);
    triangulation->getEdgeVertex(linkId, 1, neighborId1);

    // connect the two uf together
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId0) {
        uf0 = j;
        break;
      }
    }
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId1) {
        uf1 = j;
        break;
      }
    }

    seedList[uf0] = UnionFind::makeUnion(seedList[uf0], seedList[uf1]);
    seedList[uf1] = seedList[uf0];
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)seedList.size(); i++) {
    seedList[i] = seedList[i]->find();
  }

  std::vector<UnionFind *>::iterator it;
  sort(seedList.begin(), seedList.end());
  it = unique(seedList.begin(), seedList.end());
  seedList.resize(distance(seedList.begin(), it));

  return (SimplexId)seedList.size();
}
