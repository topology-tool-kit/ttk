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

namespace ttk
{
   namespace ftr
   {
      class FTRGraph : public Debug
      {
        protected:
         void *         inputData_, *outputData_;
         Triangulation *mesh_;

        public:
         FTRGraph();

         virtual ~FTRGraph();

         /// build the package.
         /// \pre If this TTK package uses ttk::Triangulation for fast mesh
         /// traversals, the function setupTriangulation() must be called on this
         /// object prior to this function, in a clearly distinct pre-processing
         /// steps. An error will be returned otherwise.
         /// \note In such a case, it is recommended to exclude
         /// setupTriangulation() from any time performance measurement.
         /// \param argment Dummy integer argument.
         /// \return Returns 0 upon success, negative values otherwise.
         template <class dataType>
         int build(const int &argument) const;

         /// Pass a pointer to an input array representing a scalarfield.
         /// The expected format for the array is the following:
         /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
         /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
         /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
         /// The array is expected to be correctly allocated.
         /// \param data Pointer to the data array.
         /// \return Returns 0 upon success, negative values otherwise.
         /// \sa setVertexNumber() and setDimensionNumber().
         inline int setInputDataPointer(void *data)
         {
            inputData_ = data;
            return 0;
         }

         /// Pass a pointer to an output array representing a scalar field.
         /// The expected format for the array is the following:
         /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
         /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
         /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
         /// The array is expected to be correctly allocated.
         /// \param data Pointer to the data array.
         /// \return Returns 0 upon success, negative values otherwise.
         /// \sa setVertexNumber() and setDimensionNumber().
         inline int setOutputDataPointer(void *data)
         {
            outputData_ = data;
            return 0;
         }

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
      };

   }  // namespace ftr
}  // namespace ttk

#include "FTRGraph_Template.h"

#endif  // FTRGRAPH_H
