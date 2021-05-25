/// \ingroup base
/// \class ttk::MultiresTriangulation
/// \author Jules Vidal
/// \date 2019.
///
/// \brief MultiresTriangulation
///
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa ProgressiveTopology
/// \sa Triangulation

#pragma once

// base code includes
#include <Geometry.h>
#include <ImplicitTriangulation.h>

#include <array>

namespace ttk {

  class MultiresTriangulation : public Debug {

  public:
    MultiresTriangulation();
    ~MultiresTriangulation();

    SimplexId getVertexNeighborAtDecimation(const SimplexId &vertexId,
                                            const int &localNeighborId,
                                            SimplexId &neighborId,
                                            int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dA(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dAB(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dABCD(const SimplexId v,
                                                  const int id,
                                                  const SimplexId shiftX,
                                                  const SimplexId shiftY,
                                                  const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dAC(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dB(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dBD(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dC(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dCD(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimation2dD(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const int decimation) const;
    SimplexId getVertexNeighborAtDecimationA(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationAB(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationABCDEFGH(const SimplexId v,
                                                    const int id,
                                                    const SimplexId shiftX,
                                                    const SimplexId shiftY,
                                                    const SimplexId shiftZ,
                                                    const int decimation) const;
    SimplexId getVertexNeighborAtDecimationABDC(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationAC(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationAE(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationAEFB(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationAEGC(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationB(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationBD(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationBF(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationBFHD(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationC(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationCD(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationCG(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationD(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationDH(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationE(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationEF(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationEFHG(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationEG(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationF(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationFH(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationG(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighborAtDecimationGH(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY,
                                              const SimplexId shiftZ,
                                              const int decimation) const;
    SimplexId getVertexNeighborAtDecimationGHDC(const SimplexId v,
                                                const int id,
                                                const SimplexId shiftX,
                                                const SimplexId shiftY,
                                                const SimplexId shiftZ,
                                                const int decimation) const;
    SimplexId getVertexNeighborAtDecimationH(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ,
                                             const int decimation) const;
    SimplexId getVertexNeighbor2dA(const SimplexId v,
                                   const int id,
                                   SimplexId shiftX,
                                   SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dB(const SimplexId v,
                                   const int id,
                                   SimplexId shiftX,
                                   SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dC(const SimplexId v,
                                   const int id,
                                   SimplexId shiftX,
                                   SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dD(const SimplexId v,
                                   const int id,
                                   SimplexId shiftX,
                                   SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dAB(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dCD(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dAC(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dBD(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY) const;
    SimplexId getVertexNeighbor2dABCD(const SimplexId v,
                                      const int id,
                                      SimplexId shiftX,
                                      SimplexId shiftY) const;

    SimplexId getVertexNeighborA(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborB(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborC(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborD(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborE(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborF(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborG(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborH(const SimplexId v,
                                 const int id,
                                 SimplexId shiftX,
                                 SimplexId shiftY,
                                 SimplexId shiftZ) const;
    SimplexId getVertexNeighborAB(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborCD(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborEF(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborGH(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborAC(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborBD(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborEG(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborFH(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborAE(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborBF(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborCG(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborDH(const SimplexId v,
                                  const int id,
                                  SimplexId shiftX,
                                  SimplexId shiftY,
                                  SimplexId shiftZ) const;
    SimplexId getVertexNeighborABDC(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborEFHG(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborAEGC(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborBFHD(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborAEFB(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborGHDC(const SimplexId v,
                                    const int id,
                                    SimplexId shiftX,
                                    SimplexId shiftY,
                                    SimplexId shiftZ) const;
    SimplexId getVertexNeighborABCDEFGH(const SimplexId v,
                                        const int id,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;

    SimplexId getInvertVertexNeighbor2dA(const SimplexId v,
                                         const SimplexId neighborId,
                                         SimplexId shiftX,
                                         SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dB(const SimplexId v,
                                         const SimplexId neighborId,
                                         SimplexId shiftX,
                                         SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dC(const SimplexId v,
                                         const SimplexId neighborId,
                                         SimplexId shiftX,
                                         SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dD(const SimplexId v,
                                         const SimplexId neighborId,
                                         SimplexId shiftX,
                                         SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dAB(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dCD(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dAC(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dBD(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY) const;
    SimplexId getInvertVertexNeighbor2dABCD(const SimplexId v,
                                            const SimplexId neighborId,
                                            SimplexId shiftX,
                                            SimplexId shiftY) const;

    SimplexId getInvertVertexNeighborA(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborB(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborC(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborD(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborE(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborF(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborG(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborH(const SimplexId v,
                                       const SimplexId neighborId,
                                       SimplexId shiftX,
                                       SimplexId shiftY,
                                       SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborAB(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborCD(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborEF(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborGH(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborAC(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborBD(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborEG(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborFH(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborAE(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborBF(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborCG(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborDH(const SimplexId v,
                                        const SimplexId neighborId,
                                        SimplexId shiftX,
                                        SimplexId shiftY,
                                        SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborABDC(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborEFHG(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborAEGC(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborBFHD(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborAEFB(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborGHDC(const SimplexId v,
                                          const SimplexId neighborId,
                                          SimplexId shiftX,
                                          SimplexId shiftY,
                                          SimplexId shiftZ) const;
    SimplexId getInvertVertexNeighborABCDEFGH(const SimplexId v,
                                              const SimplexId neighborId,
                                              SimplexId shiftX,
                                              SimplexId shiftY,
                                              SimplexId shiftZ) const;
    SimplexId
      getInvertedVertexNeighborABCDEFGH(const SimplexId v,
                                        const int id,
                                        const SimplexId shiftX,
                                        const SimplexId shiftY,
                                        const SimplexId shiftZ,
                                        SimplexId &invertedLocalNeighbor) const;
    void getImpactedVerticesA(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesB(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesC(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesD(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesE(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesF(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesG(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesH(std::array<SimplexId, 3> &p,
                              SimplexId &localNeighborId0,
                              SimplexId &localNeighborId1) const;
    void getImpactedVerticesAB(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesEF(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesCD(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesGH(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesAC(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesEG(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesAE(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesCG(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesBD(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesFH(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesBF(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesDH(std::array<SimplexId, 3> &p,
                               SimplexId &localNeighborId0,
                               SimplexId &localNeighborId1) const;
    void getImpactedVerticesABDC(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesEFHG(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesAEFB(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesGHDC(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesAEGC(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesBFHD(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesABCDEFGH(std::array<SimplexId, 3> &p,
                                     SimplexId &localNeighborId0,
                                     SimplexId &localNeighborId1) const;

    void getImpactedVertices2dABCD(std::array<SimplexId, 3> &p,
                                   SimplexId &localNeighborId0,
                                   SimplexId &localNeighborId1) const;
    void getImpactedVertices2dAB(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVertices2dCD(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVertices2dBD(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVertices2dAC(std::array<SimplexId, 3> &p,
                                 SimplexId &localNeighborId0,
                                 SimplexId &localNeighborId1) const;
    void getImpactedVerticesError(const int prev_decim,
                                  const std::array<SimplexId, 3> &p) const;

    void getInvertedLocalNeighborA(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborB(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborC(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborD(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborE(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborF(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborG(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborH(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborAB(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborEF(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborCD(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborGH(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborAC(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborEG(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborAE(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborCG(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborBD(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborFH(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborBF(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborDH(SimplexId id,
                                    SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborABDC(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborEFHG(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborAEFB(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborGHDC(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborAEGC(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;
    void getInvertedLocalNeighborBFHD(SimplexId id,
                                      SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedVertexNeighbor2dABCD(const SimplexId v,
                                      const int id,
                                      const SimplexId shiftX,
                                      const SimplexId shiftY,
                                      SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dA(SimplexId id,
                                  SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dB(SimplexId id,
                                  SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dC(SimplexId id,
                                  SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dD(SimplexId id,
                                  SimplexId &invertedLocalNeighbor) const;
    SimplexId
      getInvertedLocalNeighbor2dAB(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dAC(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dBD(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;

    SimplexId
      getInvertedLocalNeighbor2dCD(SimplexId id,
                                   SimplexId &invertedLocalNeighbor) const;

    int getInvertVertexNeighbor(const SimplexId &vertexId,
                                const SimplexId &neighborId,
                                SimplexId &localNeighborId) const;
    int getVertexNeighbor(const SimplexId &vertexId,
                          const int &localNeighborId,
                          SimplexId &neighborId) const;
    int getInteriorInvertedVertexNeighbor(SimplexId,
                                          SimplexId,
                                          SimplexId &,
                                          SimplexId &) const;
    bool areVerticesNeighbors(const SimplexId, const SimplexId) const;
    bool isBoundaryImpacted(SimplexId) const;
    SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const;
    void vertexToPosition2d(const SimplexId vertex,
                            std::array<SimplexId, 3> &p) const;
    void vertexToPosition(const SimplexId vertex,
                          std::array<SimplexId, 3> &p) const;
    SimplexId localToGlobalVertexId(const SimplexId localId) const;
    int getVertexBoundaryIndex(const SimplexId) const;

    int getDimensionality() const {
      return dimensionality_;
    }
    void setTriangulation(ImplicitTriangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_) {
        dimensionality_ = triangulation_->getDimensionality();
        std::vector<int> gridDim(3);
        triangulation_->getGridDimensions(gridDim);
        for(int i = 0; i < 3; i++) {
          gridDimensions_[i] = gridDim[i];
          nbvoxels_[i] = gridDimensions_[i] - 1;
        }

        if(dimensionality_ == 2) {
          if(gridDimensions_[0] == 1) {
            Di_ = 1;
            Di_ = 2;
          } else if(gridDimensions_[1] == 1) {
            Di_ = 0;
            Dj_ = 2;
          } else {
            Di_ = 0;
            Dj_ = 1;
          }
          vertexNumber_ = gridDimensions_[Di_] * gridDimensions_[Dj_];
          vshift_[0] = gridDimensions_[Di_];
        } else if(dimensionality_ == 1) {
          for(int k = 0; k < 3; k++) {
            if(gridDimensions_[k] > 1) {
              Di_ = k;
              break;
            }
          }
          vertexNumber_ = gridDimensions_[Di_];
        } else if(dimensionality_ == 3) {
          Di_ = 0;
          Dj_ = 1;
          Dk_ = 2;
          vertexNumber_ = gridDimensions_[Di_] * gridDimensions_[Dj_]
                          * gridDimensions_[Dk_];
          vshift_[0] = gridDimensions_[Di_];
          vshift_[1] = gridDimensions_[Dj_] * gridDimensions_[Di_];
        } else {
          printErr("Wrong dimensionality");
        }
        this->preconditionVerticesInternal();
        this->computeCoarsestDecimationLevel();

      } else {
        printErr("Empty input triangulation !");
      }
    }

    /**
     * @brief Precondition vertices
     *
     * Fill the vertexPositions_ and the vertexCoords_ arrays with the
     * VertexPosition and the 3D coordinates on the grid for every
     * vertex
     */
    int preconditionVerticesInternal();
    int computeVerticesShifts();

    inline int getDecimatedVertexNumber() const {
      return decimatedVertexNumber_;
    }
    int getVertexNumber() const {
      return vertexNumber_;
    }
    inline int getDecimation() const {
      return decimation_;
    }
    inline int getDecimationLevel() const {
      return decimationLevel_;
    }

    inline void setDecimationLevel(int decimationLevel) {
      decimationLevel_ = decimationLevel;
      decimation_ = pow(2, decimationLevel_);
      computeDecimatedDimensions();
      computeVerticesShifts();
    }

    bool isInTriangulation(const SimplexId vertexId) const;
    void
      getImpactedVertices(SimplexId vertexId, SimplexId v0[3], SimplexId v1[3]);

    void computeDecimatedDimensions() {
      int xDim = gridDimensions_[0];
      int yDim = gridDimensions_[1];
      int zDim = gridDimensions_[2];
      if(decimation_ > 1) {
        xDim = ((xDim - 1) % decimation_) ? ((xDim - 1) / decimation_ + 2)
                                          : ((xDim - 1) / decimation_ + 1);

        yDim = ((yDim - 1) % decimation_) ? ((yDim - 1) / decimation_ + 2)
                                          : ((yDim - 1) / decimation_ + 1);

        zDim = ((zDim - 1) % decimation_) ? ((zDim - 1) / decimation_ + 2)
                                          : ((zDim - 1) / decimation_ + 1);
      }
      gridDecimatedDimensions_[0] = xDim;
      gridDecimatedDimensions_[1] = yDim;
      gridDecimatedDimensions_[2] = zDim;
      decimatedVertexNumber_ = xDim * yDim * zDim;
    }

    void computeCoarsestDecimationLevel();

    int RL_to_DL(int rl);

    int DL_to_RL(int dl);

    std::vector<int> getGridDimensions() const {
      std::vector<int> dimensions(3);
      dimensions[0] = gridDimensions_[0];
      dimensions[1] = gridDimensions_[1];
      dimensions[2] = gridDimensions_[2];
      return dimensions;
    }

    ImplicitTriangulation *getTriangulation() {
      return triangulation_;
    }

    char localNeighborId(SimplexId neighborId, SimplexId vertexId);

    int setDebugLevel(const int &debugLevel) override {
      debugLevel_ = debugLevel;
      return 0;
    }
    std::vector<SimplexId> getExtendedStar(const SimplexId &vertexId) const;
    void findBoundaryRepresentatives(
      std::vector<SimplexId> &boundaryRepresentatives);

  protected:
  private:
    enum class VertexPosition : char {
      // a--------b

      LEFT_CORNER_1D, // a
      RIGHT_CORNER_1D, // b
      CENTER_1D,
      // total: 3 1D cases

      // a--------b
      // |        |
      // |        |
      // |        |
      // c--------d

      // 2D corners
      TOP_LEFT_CORNER_2D, // a
      TOP_RIGHT_CORNER_2D, // b
      BOTTOM_LEFT_CORNER_2D, // c
      BOTTOM_RIGHT_CORNER_2D, // d
      // 2D edges
      TOP_EDGE_2D, // ab
      BOTTOM_EDGE_2D, // cd
      LEFT_EDGE_2D, // ac
      RIGHT_EDGE_2D, // bd
      // 2D central strip
      CENTER_2D,
      // total: 9 2D cases

      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // 3D corners
      TOP_LEFT_FRONT_CORNER_3D, // a
      TOP_RIGHT_FRONT_CORNER_3D, // b
      BOTTOM_LEFT_FRONT_CORNER_3D, // c
      BOTTOM_RIGHT_FRONT_CORNER_3D, // d
      TOP_LEFT_BACK_CORNER_3D, // e
      TOP_RIGHT_BACK_CORNER_3D, // f
      BOTTOM_LEFT_BACK_CORNER_3D, // g
      BOTTOM_RIGHT_BACK_CORNER_3D, // h
      // 3D edges
      TOP_FRONT_EDGE_3D, // ab
      BOTTOM_FRONT_EDGE_3D, // cd
      LEFT_FRONT_EDGE_3D, // ac
      RIGHT_FRONT_EDGE_3D, // bd
      TOP_BACK_EDGE_3D, // ef
      BOTTOM_BACK_EDGE_3D, // gh
      LEFT_BACK_EDGE_3D, // eg
      RIGHT_BACK_EDGE_3D, // fh
      TOP_LEFT_EDGE_3D, // ae
      TOP_RIGHT_EDGE_3D, // bf
      BOTTOM_LEFT_EDGE_3D, // cg
      BOTTOM_RIGHT_EDGE_3D, // dh
      // 3D faces
      FRONT_FACE_3D, // abcd
      BACK_FACE_3D, // efgh
      TOP_FACE_3D, // abef
      BOTTOM_FACE_3D, // cdgh
      LEFT_FACE_3D, // aceg
      RIGHT_FACE_3D, // bdfh
      // 3D central part
      CENTER_3D,
      // total: 27 3D cases
    };

    // for every vertex, its position on the grid
    std::vector<VertexPosition> vertexPositions_{};
    // for every vertex, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> vertexCoords_{};
    // for every vertex, the corresponding shifts
    std::vector<std::array<SimplexId, 3>> vertexShifts_{};

    int dimensionality_;
    int decimation_;
    int vertexNumber_;
    int decimatedVertexNumber_;
    int decimationLevel_;
    int coarsestDL_;
    int gridDimensions_[3];
    int gridDecimatedDimensions_[3];
    int nbvoxels_[3];
    int Di_, Dj_, Dk_;
    int vshift_[2];

    ImplicitTriangulation *triangulation_;
  };
} // namespace ttk

inline bool ttk::MultiresTriangulation::areVerticesNeighbors(
  const SimplexId vertexId0, const SimplexId vertexId1) const {
  SimplexId neighborNumber = getVertexNeighborNumber(vertexId0);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId;
    getVertexNeighbor(vertexId0, i, neighborId);
    if(neighborId == vertexId1) {
      return true;
    }
  }
  return false;
}
