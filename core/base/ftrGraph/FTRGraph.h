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

#ifndef _FTRGRAPH_H
#define _FTRGRAPH_H

// base code includes
#include <Triangulation.h>

#include "DataTypes.h"
#include "DynamicGraph.h"
#include "Graph.h"
#include "Propagation.h"
#include "Scalars.h"
#include "FTRCommon.h"

namespace ttk
{
   namespace ftr
   {
      template<typename ScalarType>
      class FTRGraph : virtual public Debug, public Allocable
      {
        private:
         // Exernal fields
         Params* const              params_;
         Triangulation*             mesh_;
         Scalars<ScalarType>* const scalars_;

         const bool needDelete_;

         // Internal fields
         Graph graph_;

        public:
         FTRGraph(Params* const params, Triangulation* mesh, Scalars<ScalarType>* const scalars);
         FTRGraph();
         virtual ~FTRGraph();

         /// build the Reeb Graph
         /// \pre If this TTK package uses ttk::Triangulation for fast mesh
         /// traversals, the function setupTriangulation() must be called on this
         /// object prior to this function, in a clearly distinct pre-processing
         /// steps. An error will be returned otherwise.
         /// \note In such a case, it is recommended to exclude
         /// setupTriangulation() from any time performance measurement.
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
         inline int setupTriangulation(Triangulation *triangulation)
         {
            mesh_ = triangulation;

            if (mesh_) {
               mesh_->preprocessVertexNeighbors();
            }

            return 0;
         }

         // Accessor on the Tree
         // ---------------------

         // TODO

         // Parameters
         // ----------

         /// The nuber of threads to be used during the computation
         /// of the reeb graph
         void setThreadNumber(const idThread nb)
         {
            params_->threadNumber_ = nb;
         }

         /// Control the verbosity of the base code
         void setDebugLevel(const int lvl)
         {
            Debug::setDebugLevel(lvl);
            params_->debugLevel_ = lvl;
         }

         /// Scalar field used to compute the Reeb Graph
         void setScalars(const void* scalars)
         {
            scalars_->setScalars((ScalarType*)scalars);
         }

         /// When several points have the same scalar value,
         /// we use simulation of simplicity to distingish between
         /// them in a morse discret geometry compliant way.
         /// This is explained in the TTK report.
         /// Set the array to use here
         void setVertexSoSoffsets(idVertex* sos)
         {
            scalars_->setOffsets(sos);
         }

        protected:

         // Build functions

         /// Find the extrema from which the local propagations will start
         void leafSearch();

         /// Swipe over the dataset, using Fibonacci Heap for the
         /// propagation and Dynamic Graph to track the connected component
         /// in oreder to construct the Graph
         void swipe();

         // Print function

         void printGraph(const int lvl) const;

         void printTime(DebugTimer& timer, const string& msg, const int lvl) const;

         // Initialize functions (virtual inherit from Allocable)
         // called automatically by the build

         void alloc() override;

         void init() override;

      };

   }  // namespace ftr
}  // namespace ttk

// Implementation
#include "FTRGraph_Template.h"

#endif  // FTRGRAPH_H
