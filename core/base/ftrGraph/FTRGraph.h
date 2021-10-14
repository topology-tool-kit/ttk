/// \ingroup base
/// \class ttk::ftr::FTRGraph
/// \author charles gueunet charles.gueunet+ttk@gmail.com
/// \date 2018-01-15
///
/// \brief TTK %FTRGraph processing package.
///
/// %FTRGraph is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkFTRGraph.cpp %for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>

// local includes
#include "DynamicGraph.h"
#include "FTRAtomicVector.h"
#include "FTRCommon.h"
#include "FTRDataTypes.h"
#include "FTRPropagation.h"
#include "FTRPropagations.h"
#include "FTRScalars.h"
#include "Graph.h"
#include "Mesh.h"

// other baseCode
#include "ScalarFieldCriticalPoints.h"

#ifndef TTK_DISABLE_FTR_LAZY
#include "FTRLazy.h"
#endif

// c++ includes
#include <set>
#include <tuple>

namespace ttk {
  namespace ftr {

    struct DynGraphs {
      // one going up, one going down
      DynamicGraph<idVertex> up{}, down{};
    };

    struct Valences {
      std::vector<valence> lower{}, upper{};
    };

    struct LocalForests {
      // one for upper link, one for lower link
      LocalForest<idVertex> up{}, down{};
    };

    struct Star {
      std::vector<idEdge> lower{}, upper{};
    };

    struct Comp {
      std::set<DynGraphNode<idVertex> *> lower{}, upper{};
    };

    template <typename ScalarType, typename triangulationType>
    class FTRGraph : public Allocable {
      // Exernal fields
      Params params_{};
      Scalars<ScalarType> scalars_{};

      // Internal fields
      Graph graph_{};
      Mesh<triangulationType> mesh_{};
      Propagations propagations_{};
      DynGraphs dynGraphs_{};
      Valences valences_{};

#ifndef TTK_DISABLE_FTR_LAZY
      Lazy lazy_{};
#endif

#ifdef TTK_ENABLE_FTR_TASK_STATS
      // Stats
      Timer sweepStart_{};
      std::vector<float> propTimes_{};
      idVertex nbProp_{};
#endif

    public:
      explicit FTRGraph(triangulationType *mesh);
      FTRGraph();

      /// build the Reeb Graph
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh
      /// traversals, the function preconditionTriangulation() must be called on
      /// this object prior to this function, in a clearly distinct
      /// pre-processing steps. An error will be returned otherwise. \note In
      /// such a case, it is recommended to exclude preconditionTriangulation()
      /// from any time performance measurement.
      void build();

      // General documentation info:
      //
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
      /// see the setInput*() functions of ttk::Triangulation). See vtkFTRGraph
      /// for further examples.
      ///
      /// \param triangulation Pointer to a valid triangulation.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa ttk::Triangulation
      //
      //
      // Developer info:
      // ttk::Triangulation is a generic triangulation representation that
      // enables fast mesh traversal, either on explicit triangulations (i.e.
      // tet-meshes) or implicit triangulations (i.e. low-memory footprint
      // implicit triangulations obtained from regular grids).
      //
      // Not all TTK packages need such mesh traversal features. If your
      // TTK package needs any mesh traversal procedure, we recommend to use
      // ttk::Triangulation as described here.
      //
      // Each call to a traversal procedure of ttk::Triangulation
      // must satisfy some pre-condition (see ttk::Triangulation for more
      // details). Such pre-condition functions are typically called from this
      // function.
      inline int preconditionTriangulation(triangulationType *triangulation) {
        mesh_.setTriangulation(triangulation);
        if(triangulation) {
          mesh_.preprocess();
        }

        return 0;
      }

      // Accessor on the graph
      // ---------------------

      Graph &&extractOutputGraph(void) {
        return std::move(graph_);
      }

      // Parameters
      // ----------

      /// The nuber of threads to be used during the computation
      /// of the reeb graph
      int setThreadNumber(const int nb) override {
        params_.threadNumber = nb;
        // Security, but do not rely on this one
        threadNumber_ = nb;
        return 0;
      }

      /// Control the verbosity of the base code
      virtual int setDebugLevel(const int &lvl) override {
        params_.debugLevel = lvl;
        return Debug::setDebugLevel(lvl);
      }

      void setParams(const Params &p) {
        params_ = p;
        threadNumber_ = params_.threadNumber;
      }

      /// Scalar field used to compute the Reeb Graph
      void setScalars(const void *scalars) {
        scalars_.setScalars((ScalarType *)scalars);
      }

      /// When several points have the same scalar value,
      /// we use simulation of simplicity to distingish between
      /// them in a morse discret geometry compliant way.
      /// This is explained in the TTK report.
      /// Set the array to use here
      void setVertexSoSoffsets(SimplexId *sos) {
        scalars_.setOffsets(sos);
      }

      DynamicGraph<idVertex> &dynGraph(const Propagation *const lp) {
        if(lp->goUp()) {
          return dynGraphs_.up;
        } else {
          return dynGraphs_.down;
        }
      }

      DynamicGraph<idVertex> &dynGraph(const bool goUp) {
        if(goUp) {
          return dynGraphs_.up;
        } else {
          return dynGraphs_.down;
        }
      }

      DynamicGraph<idVertex> &dynGraphOpposite(const Propagation *const lp) {
        if(lp->goUp()) {
          return dynGraphs_.down;
        } else {
          return dynGraphs_.up;
        }
      }

    protected:
      // Build functions

      // classify critical points, marks saddle in join/split vectors
      // and add min/max or both as leaves.
      void criticalSearch();

      /// Launch the sweep algorithm, but adapted to growth
      /// locally from each seed (min and/or max)
      /// See: grwothFromSeed.
      void sweepFrowSeeds();

      // launch a sequential sweep on the whole mesh.
      void sweepSequential();

      // Print function (FTRGraphPrint)

      std::string printMesh(void) const;

      std::string printEdge(const idEdge edgeId,
                            const Propagation *const localProp) const;

      std::string printTriangle(const idCell cellId,
                                const Propagation *const localProp) const;

      void printGraph(const int verbosity) const;

      void printTime(Timer &timer, const std::string &msg) const;

      // Initialize functions (virtual inherit from Allocable)
      // called automatically by the build

      void alloc() override;

      void init() override;

    private:
      /// Local propagation for the vertex seed, using BFS with a priority queue
      /// localProp.
      /// This will process all the area corresponding to one connected
      /// component of level set. When a 2 Saddle is met, this function wait it
      /// has been completely visited and then continue. When a 1 Saddle is met,
      /// we split the local propagation with a BFS to continue locally. if arc
      /// is supplied, this arc will be used for the growth NOTE: use an
      /// insertion/deletion list to add lazyness on DynGraph
      void growthFromSeed(const idVertex seed,
                          Propagation *localProp,
                          idSuperArc currentArc = nullSuperArc);

      // growth like in the original algorithm by maintaining one dynamic graph.
      // the direction of the growth: increasing scalar value.
      void growthSequential(const idVertex begin, const idVertex stop);

      /// visit the star of v and fill the two vectors in Star,
      void visitStar(const Propagation *const localProp, Star &star) const;

      /// Consider edges ending at the vertex v, one by one,
      /// and find their corresponding components in the current
      /// preimage graph, each representing a component.
      /// \ret the set of uniques representing components
      std::set<DynGraphNode<idVertex> *>
        lowerComps(const std::vector<idEdge> &finishingEdges,
                   const Propagation *const localProp);

      /// Symetric to lowerComps
      /// \ref lowerComps
      std::set<DynGraphNode<idVertex> *>
        upperComps(const std::vector<idEdge> &startingEdges,
                   const Propagation *const localProp);

      bool checkStop(const std::vector<DynGraphNode<idVertex> *> &lowerComp);

      // visit these edges neighborhood to understand the local connectivity
      // return the number of component in lower/upper link
      std::pair<valence, valence> getLinkNbCC(const idVertex curVert,
                                              LocalForests &localForests,
                                              const VertCompFN &comp);

      /// update (locally) the preimage graph (dynGraph) from that
      /// of immediately before f(v) to that of immediately after f(v).
      /// (v is localGrowth->getCurVertex())
      // this update is made using the arc: curArc
      void updatePreimage(const Propagation *const localProp,
                          const idSuperArc curArc);

      /// update the dynamicGraph by adding (if needed) a new edge corresponding
      /// to the starting cell cellId
      void updatePreimageStartCell(const orderedTriangle &oTriangle,
                                   const Propagation *const localProp,
                                   const idSuperArc curArc);

      /// update the dynamicGraph by moving (if needed) the corresponding to the
      /// current visited cell cellId
      void updatePreimageMiddleCell(const orderedTriangle &oTriangle,
                                    const Propagation *const localProp,
                                    const idSuperArc curArc);

      /// update the dynamicGraph by removing (if needed) the edge corresponding
      /// to the last visit of the cell cellId
      void updatePreimageEndCell(const orderedTriangle &oTriangle,
                                 const Propagation *const localProp,
                                 const idSuperArc curArc);

#ifndef TTK_DISABLE_FTR_LAZY

      // like updatePreimage but only impacte the lazy list in propagation
      // and not the global DG
      void lazyUpdatePreimage(Propagation *const localProp,
                              const idSuperArc curArc);

      // updatePreimageStart but in lazy, impact lazy lists and not the DG
      void updateLazyStart(const orderedTriangle &oTriangle,
                           Propagation *const localProp,
                           const idSuperArc curArc);

      // updatePreimageMiddle but in lazy, impact lazy lists and not the DG
      void updateLazyMiddle(const orderedTriangle &oTriangle,
                            Propagation *const localProp,
                            const idSuperArc curArc);

      // updatePreimageEnd but in lazy, impact lazy lists and not the DG
      void updateLazyEnd(const orderedTriangle &oTriangle,
                         Propagation *const localProp,
                         const idSuperArc curArc);

      // impacte the global DG
      void updateLazyAdd(const Propagation *const localProp,
                         const linkEdge edge,
                         const idSuperArc arc);

      void updateLazyDel(const Propagation *const localProp,
                         const linkEdge edge,
                         const idSuperArc arc);

      // aply modifications from lazy lists to the global DG (for arc a)
      void lazyApply(Propagation *const locapProp, const idSuperArc a);

#endif

      /// update the current arc of the dynGraph subtree of seed with curArc (on
      /// the component going through neigEdge)
      void updateDynGraphCurArc(const idVertex seed,
                                const idEdge neigEdge,
                                const idSuperArc curArc,
                                const Propagation *const localProp);

      /// update the current arc of the dynGraph subtree of seed with curArc
      /// (all edges crossing the level set)
      void updateDynGraphCurArc(const idVertex seed,
                                const idSuperArc curArc,
                                const Propagation *const localProp);

      /// update the skeleton structure
      /// \ret the nodeId of the current saddle/max
      idNode updateReebGraph(const idSuperArc currentArc,
                             const Propagation *const localProp);

      /// local growth replacing the global sort
      /// Add vertices above the current one in the propagation,
      /// return true if vertices aboves the current one have been found.
      /// Note, these vertices may not have been added if already marked as in
      /// the propagation.
      void localGrowth(Propagation *const localProp,
                       const std::vector<idEdge> &upperEdges);

      // Check if the current vertex which is on a Join saddle come from the
      // last growth touching this saddle
      bool checkLast(Propagation *const localProp,
                     const std::vector<idEdge> &lowerStarEdges);

      // At a join saddle, merge local propagations coming here
      // and close remiang opened arcs.
      // Remove duplicate on the saddleVertex (only)
      // return the number of visible arcs merging
      idSuperArc
        mergeAtSaddle(const idNode saddleId,
                      Propagation *localProp,
                      const std::set<DynGraphNode<idVertex> *> &lowerComp);

      // At a join saddle, close onped arcs only
      // do not touch local propagations
      // return the number of visible arcs merging
      idSuperArc
        mergeAtSaddle(const idNode saddleId,
                      const std::set<DynGraphNode<idVertex> *> &lowerComp);

      // At a split saddle, assign new arcs at each CC in the DynGraph,
      // and launch a new propagation taking care of these arcs simultaneously
      // if hidden is true, new arcs are created hidden
      void splitAtSaddle(Propagation *const localProp,
                         const std::set<DynGraphNode<idVertex> *> &upperComp,
                         const bool hidden = false);

      // Retrun one triangle by upper CC of the vertex v
      std::set<idCell> upCCtriangleSeeds(const idVertex v,
                                         const Propagation *const localProp);

      // bfs on triangles/edges in the neighborhood of the saddle to mark
      // each cc with a distinct identifier.
      void bfsSeed(const std::size_t idt,
                   const valence idcc,
                   std::vector<idCell> &triangles,
                   std::vector<valence> &cc,
                   const Propagation *const localProp);

      // bfs on triangles/edges crossing the level set at saddle, starting
      // at seed. upper vertices encountred are added to newLocalProp
      // : saddle is the starting saddle,
      // : seed is at first call the first triangle of this bfs (will change
      // during recursive call)
      // : propagation is the propagation to fill with upper vertcices of
      // visited edges i
      // : arc is the arc id used to marks the segmentation of lower
      // vertices of visited edges (it has newLocalProp associated
      void bfsPropagation(const idVertex saddle,
                          const idCell seed,
                          Propagation *const newLocalProp,
                          const idSuperArc arc);

      // visit a vertex in terms of segmantation and history,
      // also check if the current arc is merging through an opposite one.
      idSuperArc visit(Propagation *const localProp, const idSuperArc curArc);

      // Tools

      // Create a new propagation starting at leaf
      Propagation *newPropagation(const idVertex leaf, const bool fromMax);

      // Compute the wieght of the edge in the dyngraph between e1 and e2.
      // This weight is the min value of the two endpoints, we use the mirror
      // array (int)
      idVertex getWeight(const orderedEdge &e1,
                         const orderedEdge &e2,
                         const Propagation *const localProp);

      // DEPRECATED
      // TODO move in Mesh if not already done

      /// On a triangle, recover the position of the current vertex to classify
      /// the triangle
      vertPosInTriangle
        getVertPosInTriangle(const orderedTriangle &oTriangle,
                             const Propagation *const localProp) const;

      // get the higher vertex: get<1>(get<1>(oTriangle)
      idVertex getEndVertexInTriangle(const orderedTriangle &oTriangle,
                                      const Propagation *const localProp) const;

      // get edge in orderer triangle between v0 and v1
      idEdge getEdgeFromOTri(const orderedTriangle oTri,
                             const idVertex v0,
                             const idVertex v1);
    };
  } // namespace ftr
} // namespace ttk

// Implementation
#include "FTRGraphPrint_Template.h"
#include "FTRGraphPrivate_Template.h"
#include "FTRGraph_Template.h"
