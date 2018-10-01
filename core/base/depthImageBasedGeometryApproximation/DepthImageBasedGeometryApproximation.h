/// \ingroup base
/// \class ttk::DepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.7.2018
///
/// \brief TTK %depthImageBasedGeometryApproximation processing package.
///
/// %DepthImageBasedGeometryApproximation is a TTK processing package that derives geomerty based on an input depth image and its corresponding camera parameters.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include                 <Wrapper.h>

using namespace std;

namespace ttk{

  class DepthImageBasedGeometryApproximation : public Debug{

    public:

      DepthImageBasedGeometryApproximation(){
          inputData_ = NULL;
      };
      ~DepthImageBasedGeometryApproximation(){};

      /// Pass a pointer to an input array containing the depth values
      inline int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }

      /// Execute the Geometry Approximation.
      template <class dataType>
        int execute(
            double* camPos,
            double* camDir,
            double* camUp,
            double* camRes,
            double* camNearFar,
            double* camHeight,

            int downsampling,

            vector<tuple<float,float,float>>& vertices,
            vector<tuple<int,int,int>>& triangles,
            vector<float>& triangleDistortions
        ) const;

    protected:
      void  *inputData_;
  };

}

// if the package is a pure template class, uncomment the following line
// #include                  <DepthImageBasedGeometryApproximation.cpp>
template <class dataType> int ttk::DepthImageBasedGeometryApproximation::execute(
    double* camPos,
    double* camDir,
    double* camUp,
    double* camRes,
    double* camNearFar,
    double* camHeight,

    int downsampling,

    vector<tuple<float,float,float>>& vertices,
    vector<tuple<int,int,int>>& triangles,
    vector<float>& triangleDistortions
) const{

    Timer t;
    int numberInitialVertices = vertices.size();
    int numberInitialTriangles = triangles.size();
    size_t step = downsampling + 1;

    size_t camResST[2] = {(size_t) camRes[0], (size_t) camRes[1]};

    dataType* depthValues = (dataType*) this->inputData_;

    // camRight = camDir x CamUp
    double camRight[3] = {
        camDir[1]*camUp[2] - camDir[2]*camUp[1],
        camDir[2]*camUp[0] - camDir[0]*camUp[2],
        camDir[0]*camUp[1] - camDir[1]*camUp[0]
    };

    double temp = sqrt( camRight[0]*camRight[0] + camRight[1]*camRight[1] + camRight[2]*camRight[2] );
    camRight[0]/=temp;
    camRight[1]/=temp;
    camRight[2]/=temp;

    double camUpTrue[3] = {
        camDir[1]*(-camRight[2]) - camDir[2]*(-camRight[1]),
        camDir[2]*(-camRight[0]) - camDir[0]*(-camRight[2]),
        camDir[0]*(-camRight[1]) - camDir[1]*(-camRight[0])
    };
    temp = sqrt( camUpTrue[0]*camUpTrue[0] + camUpTrue[1]*camUpTrue[1] + camUpTrue[2]*camUpTrue[2] );
    camUpTrue[0]/=temp;
    camUpTrue[1]/=temp;
    camUpTrue[2]/=temp;

    // cout<<camUpTrue[0]<<" "<<camUpTrue[1]<<" "<<camUpTrue[2]<<endl;
    // cout<<camRight[0]<<" "<<camRight[1]<<" "<<camRight[2]<<endl;
    // cout<<camDir[0]<<" "<<camDir[1]<<" "<<camDir[2]<<endl;

    double camSize[2] = {
        camRes[0]/camRes[1]*camHeight[0],
        camHeight[0]
    };

    // Compute Index Map
    size_t n = camResST[0] * camResST[1];

    vector<int> pixelIndexVertexIndexMap;
    pixelIndexVertexIndexMap.resize( n );
    size_t numberNewVertices = 0;
    {
        for(size_t i=0; i<n; i++)
            pixelIndexVertexIndexMap[i] = depthValues[i] > 0.99 ? -1 : numberNewVertices++;
    }

    // -------------------------------------------------------------------------
    // Create Vertices
    // -------------------------------------------------------------------------
    {
        double pixelWidthWorld = camSize[0];
        double pixelHeightWorld = camSize[1];

        // make room for new vertices
        vertices.resize( numberInitialVertices + numberNewVertices );

        double delta = camNearFar[1]-camNearFar[0];

        // compute vertex positions
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t y=0; y<camResST[1]; y+=step){
            size_t yD = y*camResST[0];
            double v = (((double)y)/camRes[1] - 0.5) * pixelHeightWorld;

            for(size_t x=0; x<camResST[0]; x+=step){
                size_t pixelIndex = x + yD;

                int vertexIndex = pixelIndexVertexIndexMap[ pixelIndex ];
                if(vertexIndex < 0) continue;

                double d = (double)(depthValues[ pixelIndex ])*delta+camNearFar[0];
                double u = ( ((double)x)/camRes[0] - 0.5) * pixelWidthWorld;
                auto& vertex = vertices[vertexIndex];

                get<0>(vertex) = camPos[0] + u*camRight[0] + v*camUpTrue[0] + d*camDir[0];
                get<1>(vertex) = camPos[1] + u*camRight[1] + v*camUpTrue[1] + d*camDir[1];
                get<2>(vertex) = camPos[2] + u*camRight[2] + v*camUpTrue[2] + d*camDir[2];

                // cout<<get<0>(vertex)<<" "<<get<1>(vertex)<<" "<<get<2>(vertex)<<endl;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Create Triangles
    // -------------------------------------------------------------------------
    {
        auto absDiff = [](dataType a, dataType b){
            return a>b ? a-b : b-a;
        };

        /* Index Structure:
        0 - 1
        |   |
        2 - 3
        */
        size_t xl = camResST[0]-step;
        size_t yl = camResST[1]-step;
        size_t yD = step*camResST[0];

        for(size_t y=0; y<yl; y+=step){
            for(size_t x=0; x<xl; x+=step){
                size_t i0 = x  + y*camResST[0];
                size_t i1 = i0 + step;
                size_t i2 = i0 + yD;
                size_t i3 = i2 + step;

                int i0Index = pixelIndexVertexIndexMap[i0];
                int i1Index = pixelIndexVertexIndexMap[i1];
                int i2Index = pixelIndexVertexIndexMap[i2];
                int i3Index = pixelIndexVertexIndexMap[i3];

                dataType i0Depth = depthValues[i0];
                dataType i1Depth = depthValues[i1];
                dataType i2Depth = depthValues[i2];
                dataType i3Depth = depthValues[i3];

                if(pixelIndexVertexIndexMap[i1]>=0 && pixelIndexVertexIndexMap[i2]>=0){
                    // Check first triangle
                    if(pixelIndexVertexIndexMap[i0]>=0){
                        triangles.push_back( make_tuple( i0Index, i1Index, i2Index ) );

                        dataType distortion = max(
                            absDiff(i0Depth,i1Depth),
                            max( absDiff(i1Depth,i2Depth), absDiff(i0Depth,i2Depth) )
                        );

                        triangleDistortions.push_back( distortion );
                    }

                    // Check second triangle
                    if(pixelIndexVertexIndexMap[i3]>=0){
                        triangles.push_back( make_tuple( i1Index, i3Index, i2Index ) );

                        dataType distortion = max(
                            absDiff(i3Depth,i1Depth),
                            max( absDiff(i1Depth,i2Depth), absDiff(i3Depth,i2Depth) )
                        );

                        triangleDistortions.push_back( distortion );
                    }
                }
            }
        }
    }

    {
        std::stringstream msg;
        msg << "[ttkDepthImageBasedGeometryApproximation] Depth Image ("<<camResST[0]<<"x"<<camResST[1]<<":"<<step<<") processed in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << std::endl;
        msg << "[ttkDepthImageBasedGeometryApproximation] " << "generated (" << (vertices.size()-numberInitialVertices) << " vertices) and (" << (triangles.size()-numberInitialTriangles) << " triangles)" << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
    }

  return 0;
}
