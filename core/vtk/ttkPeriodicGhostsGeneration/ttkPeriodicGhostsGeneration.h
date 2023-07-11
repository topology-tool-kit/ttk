/// \ingroup vtk
/// \class ttkPeriodicGhostsGeneration
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date June 2023.
///
/// \brief TTK VTK-filter that generates an outside ghost layer for periodic
/// implicit
///  grids when used in a distributed context.
///
/// \param Input Scalar field on an vtkImageData.
/// \param Output Scalar field on an vtkImageData.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::PeriodicGhostsGeneration
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkPeriodicGhostsGenerationModule.h>

// VTK Includes
#include <RegularGridTriangulation.h>
#include <Triangulation.h>
#include <ttkAlgorithm.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStructuredPoints.h>

#ifdef TTK_ENABLE_MPI
#include <cstring>
#include <vtkCommunicator.h>
#include <vtkExtentTranslator.h>
#include <vtkExtractVOI.h>
#endif

namespace ttk {

  namespace periodicGhosts {
    struct partialGlobalBound {
      unsigned char isBound{0};
      double x{0};
      double y{0};
      double z{0};
    };
  } // namespace periodicGhosts
} // namespace ttk

class TTKPERIODICGHOSTSGENERATION_EXPORT ttkPeriodicGhostsGeneration
  : public ttkAlgorithm {

private:
#ifdef TTK_ENABLE_MPI
  std::array<int, 6> outExtent_; // Global extent size before computation
  std::array<int, 6> inExtent_; // Global extent size after computation
  std::array<double, 6> boundsWithoutGhosts_; // Local Bounds without ghosts
  std::array<double, 6> globalBounds_;
  std::array<double, 3> origin_;
  std::array<double, 3> spacing_;
  bool isOutputExtentComputed_{false};
  std::array<ttk::periodicGhosts::partialGlobalBound, 6> localGlobalBounds_;
  std::vector<int> neighbors_;
  std::vector<std::array<ttk::SimplexId, 6>> neighborVertexBBoxes_;
  std::array<unsigned char, 6> isBoundaryPeriodic_{};
#endif
public:
  static ttkPeriodicGhostsGeneration *New();
  vtkTypeMacro(ttkPeriodicGhostsGeneration, ttkAlgorithm);

  ttkPeriodicGhostsGeneration();
  ~ttkPeriodicGhostsGeneration() override = default;

#ifdef TTK_ENABLE_MPI

  int RequestUpdateExtent(vtkInformation *ttkNotUsed(request),
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector) override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVectors,
                         vtkInformationVector *outputVector) override;
  /**
   * This method is called in GetTriangulation, if the triangulation is
   * periodic, to create ghosts specific to dealing with this type of
   * triangulation. This may add points to the dataset of a process and
   * therefore invalidates the triangulation object taken as parameter here.
   */
  int MPIPeriodicGhostPipelinePreconditioning(vtkImageData *imageIn,
                                              vtkImageData *imageOut);

  inline std::vector<int> getNeighbors() {
    std::sort(neighbors_.begin(), neighbors_.end());
    return neighbors_;
  }

  inline std::array<unsigned char, 6> getIsBoundaryPeriodic() {
    return isBoundaryPeriodic_;
  }

#endif

protected:
  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

#ifdef TTK_ENABLE_MPI

  /**
   * @brief  Computes several pieces of information necessary for predicting the
   * global output extent.
   *
   * @return int Returns 1 if success
   */
  int ComputeOutputExtent();

  /**
   * @brief Converts extracted sub-extents of a vtkImageData to vtkCharArrays
   * and exchanges them between processes using MPI_SendRecv.
   *
   * @tparam matchesSize
   * @tparam metaDataSize
   * @param imageIn : vtkImageData input data set
   * @param charArrayBoundaries: array that will store the vtkCharArrays
   * @param charArrayBoundariesMetaData : array that will store the meta data of
   * the vtkCharArrays (such as the boundaries and processes involved in the
   * exchange)
   * @param matches : vector of boundary matches (1D, 2D or 3D)
   * @param charArrayBoundariesReceived : array that will store the received
   * vtkCharArrays
   * @param charArrayBoundariesMetaDataReceived : array that will store the meta
   * data of the received vtkCharArrays
   * @param dim Dimension of the boundaries (1, 2 or 3)
   * @return int Returns 1 upon success
   */
  template <int matchesSize, int metaDataSize>
  int MarshalAndSendRecv(
    vtkImageData *imageIn,
    std::vector<std::vector<vtkSmartPointer<vtkCharArray>>>
      &charArrayBoundaries,
    std::vector<std::vector<std::array<ttk::SimplexId, metaDataSize>>>
      &charArrayBoundariesMetaData,
    std::vector<std::array<ttk::SimplexId, matchesSize>> &matches,
    std::vector<vtkSmartPointer<vtkCharArray>> &charArrayBoundariesReceived,
    std::vector<std::array<ttk::SimplexId, metaDataSize>>
      &charArrayBoundariesMetaDataReceived,
    int dim);

  /**
   * @brief Merges the data arrays of a data set and a slice of a data set.
   *
   * A slice is a vtkImageData for which one of the dimensions is equal to 1.
   *
   * @param image Input vtkImageData
   * @param slice Slice input vtkImageData
   * @param mergedImage Output vtkImageData
   * @param direction Direction along which the merge should be done
   * @return int Returns 1 upon success
   */
  int MergeImageAndSlice(vtkImageData *image,
                         vtkImageData *slice,
                         vtkImageData *mergedImage,
                         int direction);

  /**
   * @brief Searches in the meta data array if a boundary matching direction as
   * been received. If so, converts the vtkCharArray to its original form
   * (vtkImageData) and merges the obtained vtkImageData with the mergedImage.
   *
   * @tparam boundaryType
   * @param metaDataReceived : meta data of the vtkCharArrays
   * @param boundariesReceived : vtkCharArrays to convert
   * @param direction : Directions characterizing the boundaries
   * @param mergeDirection : direction along which to merge the data sets
   * @param mergedImage: merged image
   * @return int Returns 1 upon success
   */
  template <typename boundaryType>
  int UnMarshalAndMerge(
    std::vector<boundaryType> &metaDataReceived,
    std::vector<vtkSmartPointer<vtkCharArray>> &boundariesReceived,
    boundaryType direction,
    int mergeDirection,
    vtkImageData *mergedImage);

  /**
   * @brief Searches in the meta data array if a boundary matching direction as
   * been received. If so, converts the vtkCharArray to its original form
   * (vtkImageData) and copies the obtained vtkImageData into mergedImage.
   *
   * @tparam boundaryType
   * @param metaDataReceived meta data of the vtkCharArray
   * @param boundariesReceived  vtkCharArrays to convert
   * @param direction Directions characterizing the boundaries
   * @param mergedImage merged image
   * @return int Returns 1 upon success
   */
  template <typename boundaryType>
  int UnMarshalAndCopy(
    std::vector<boundaryType> &metaDataReceived,
    std::vector<vtkSmartPointer<vtkCharArray>> &boundariesReceived,
    boundaryType direction,
    vtkImageData *mergedImage);

  /**
  * @brief Merges two data arrays along the given direction.

  * The data array from the slice has one dimension equal to 1.
  *
  * @param imageArray : image array
  * @param sliceArray : slice array
  * @param currentArray : output array
  * @param direction : direction along which to merge
  * @param dims : dimension of the vtkImageData image
  * @param ghostValue : value to give the ghost simplices
  * @param numberOfSimplices : number of simplices in currentArray
  * @param numberOfTuples : number of tuple in the sliceArray
  * @return int Returns 1 upon success
  */
  int MergeDataArrays(vtkDataArray *imageArray,
                      vtkDataArray *sliceArray,
                      vtkSmartPointer<vtkDataArray> &currentArray,
                      int direction,
                      int dims[3],
                      unsigned char ghostValue,
                      ttk::SimplexId numberOfSimplices,
                      ttk::SimplexId numberOfTuples);

#endif

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
