/// \ingroup base
/// \class ttk::ImageMetrics
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.10.2018
///
/// \brief TTK %imageMetrics processing package.
///
/// %ImageMetrics is a TTK processing package that computes the ADD metric for two images.

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk{
    class ImageMetrics : public Debug{
        public:
            ImageMetrics(){};
            ~ImageMetrics(){};

            // Computes ADD metric
            template <class dataType> int ADD(
                // Input
                dataType* pixelsA,
                dataType* pixelsB,
                size_t n,

                // Output
                double& add
            ) const;
    };
}

template <class dataType> int ttk::ImageMetrics::ADD(
    // Input
    dataType* pixelsA,
    dataType* pixelsB,
    size_t n,

    // Output
    double& add
) const {
    add = 0;

    for(size_t i=0; i<n; i++){
        add += abs( pixelsA[i] - pixelsB[i] );
    }
    add /= n;

    return 0;
};