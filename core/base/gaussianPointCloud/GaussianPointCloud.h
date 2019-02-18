/// \ingroup base
/// \class ttk::GaussianPointCloud
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK %gaussianPointCloud processing package that generates an Icosphere.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Wrapper.h>
#include <map>

using namespace std;

namespace ttk{

    class GaussianPointCloud : public Debug{

        public:

            GaussianPointCloud();
            ~GaussianPointCloud();

            // Execute the geometry approximation.
            int generate(
                // Input
                size_t subdivisions,
                float radius,
                float* center,

                // Output
                vector< tuple<float,float,float> >& vertices,
                vector< tuple<long long,long long,long long> >& triangles
            ) const;
    };
}