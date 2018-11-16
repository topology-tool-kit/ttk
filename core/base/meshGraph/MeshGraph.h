/// \ingroup base
/// \class ttk::MeshGraph
/// \author Wiebke Koepp (wiebke.koepp@gmail.com) and Jonas Lukasczyk (jl@jluk.de)
/// \date 1.11.2018
///
/// \brief TTK %meshGraph processing package.
///
/// %MeshGraph is a TTK processing package that TODO.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk{

    class MeshGraph : public Debug{

        public:

            MeshGraph(){};
            ~MeshGraph(){};

            // Execute the geometry approximation.
            template <class dataType> int execute(
                // Input
                const float* inputVertices,
                const dataType* inputTopology,
                const double* inputVertexWidths,
                size_t nVertices,
                size_t nEdges,
                size_t nSubdivisions,

                // Output
                float* outputVertices,
                dataType* outputTopology
            ) const;
    };
}

template <class dataType> int ttk::MeshGraph::execute(
    // Input
    const float* inputVertices,
    const dataType* inputTopology,
    const double* inputVertexWidths,
    size_t nVertices,
    size_t nEdges,
    size_t nSubdivisions,

    // Output
    float* outputVertices,
    dataType* outputTopology
) const{

    Timer t;

    // Print Input
    {
        stringstream msg;
        msg << "[ttkMeshGraph] Computing Mesh for graph with" << endl
            << "[ttkMeshGraph]  - "<< nVertices << " vertices" << endl
            << "[ttkMeshGraph]  - "<< nEdges << " edges" << endl;
        dMsg(cout, msg.str(), timeMsg);
    }



    // =========================================================================
    // Comment Level 0
    // =========================================================================

    // -------------------------------------------------------------------------
    // Comment Level 1
    // -------------------------------------------------------------------------

    // Comment Level 2

    // Print performance
    {
        stringstream msg;
        msg << "[ttkMeshGraph] Layout for  computed in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 0;
}
