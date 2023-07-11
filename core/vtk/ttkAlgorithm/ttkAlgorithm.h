/// \defgroup vtk vtk
/// \brief The Topology ToolKit - VTK wrapping code for the processing
/// packages.
/// @{
/// \class ttkAlgorithm
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2019.
///
/// \brief Baseclass of all VTK filters that wrap ttk modules.
///
/// This is an abstract vtkAlgorithm that provides standardized input/output
/// management for VTK wrappers of ttk filters. The class also provides a static
/// method to retrieve a ttk::Triangulation of a vtkDataSet.

#pragma once

// VTK Module
#include <ttkAlgorithmModule.h>

// VTK Includes
#include <vtkAlgorithm.h>
class vtkDataSet;
class vtkInformation;
class vtkInformationIntegerKey;

// Base Includes
#include <Debug.h>

namespace ttk {
  class Triangulation;
}

class TTKALGORITHM_EXPORT ttkAlgorithm : public vtkAlgorithm,
                                         virtual public ttk::Debug {
private:
  int ThreadNumber{1};
  bool UseAllCores{true};

public:
  static ttkAlgorithm *New();
  vtkTypeMacro(ttkAlgorithm, vtkAlgorithm);

  /**
   * Updates the number of threads of the base class based on the settings of
   * the VTK wrapper.
   */
  void UpdateThreadNumber() {
    // update ttk::BaseClass member
    this->setThreadNumber(this->UseAllCores ? ttk::OsCall::getNumberOfCores()
                                            : this->ThreadNumber);
    this->Modified();
  }

  /**
   * Explicitly sets the maximum number of threads for the base code
   * (overridden by UseAllCores member).
   */
  void SetThreadNumber(int threadNumber) {
    this->ThreadNumber = threadNumber;
    this->UpdateThreadNumber();
  }

  /**
   * Controls if the base code should use all available cores.
   */
  void SetUseAllCores(bool useAllCores) {
    this->UseAllCores = useAllCores;
    this->UpdateThreadNumber();
  }

  /**
   * Controls the debug level used by algorithms that are invoked by the VTK
   * wrapper.
   */
  void SetDebugLevel(int debugLevel) {
    this->setDebugLevel(debugLevel); // from ttk::Debug
    this->Modified();
  }

  /**
   * Set the cache size of the compact triangulation.
   */
  void SetCompactTriangulationCacheSize(float cacheSize) {
    this->CompactTriangulationCacheSize = cacheSize;
    this->Modified();
  }

  /// This method retrieves an optional array to process.
  /// The logic of this method is as follows:
  ///   - if \p enforceArrayIndex is set to true, this method will try to
  ///   retrieve the optional array to process by its identifier \p arrayIndex
  ///   (see GetInputArrayToProcess() and SetInputArrayToProcess())
  ///   - if \p enforceArrayIndex is set to false, this method will try to
  ///   retrieve the optional array to process by its name \p arrayName.
  ///
  /// In both cases, this information will be retrieved on port \p inputPort.
  ///
  vtkDataArray *GetOptionalArray(const bool &enforceArrayIndex,
                                 const int &arrayIndex,
                                 const std::string &arrayName,
                                 vtkDataSet *const inputData,
                                 const int &inputPort = 0);

  /**
   * Returns a string containing the name of the corresponding offset
   * field from a given scalar field
   */
  static std::string GetOrderArrayName(vtkDataArray *const array);

  /**
   * Retrieves an offset field from the given scalar field \p sfArray
   * or generates one, either disambiguated with the implicit vertex
   * identifier field, or with a user-provided offset field through
   * the \p enforceArrayIndex parameter and the \p arrayIndex. The
   * generated sorted offset field is then attached to the input
   * vtkDataset \p inputData.
   */
  vtkDataArray *GetOrderArray(vtkDataSet *const inputData,
                              const int scalarArrayIdx,
                              const int orderArrayIdx = 0,
                              const bool enforceOrderArrayIdx = false);

  /**
   * Retrieve an identifier field and provides a ttk::SimplexId
   * pointer to the underlying buffer.
   *
   * Use the same parameters as GetOptionalArray to fetch the VTK data
   * array.
   *
   * Fills the vector \p spareStorage if the VTK data array is not a
   * ttkSimplexIdTypeArray. This vector should have a lifetime of a
   * least the filter's RequestData method.
   */
  ttk::SimplexId *
    GetIdentifierArrayPtr(const bool &enforceArrayIndex,
                          const int &arrayIndex,
                          const std::string &arrayName,
                          vtkDataSet *const inputData,
                          std::vector<ttk::SimplexId> &spareStorage,
                          const int inputPort = 0,
                          const bool printErr = true);

  /**
   * This method retrieves the ttk::Triangulation of a vtkDataSet.
   *
   * Note, this method initializes a triangulation if one does not exist
   * already, and updates the triangulation if the connectivity of vtkDataSet
   * changed since the last retrieval.
   *
   * In the current implementation, a pointer to the triangulation of a
   * vtkDataSet is stored as a special field data array (ttkTriangulationArray)
   * of the vtkDataSet. First, the method checks if such an array exists. If
   * this is not the case, then the method initializes a triangulation for the
   * given vtkDataSet and adds a corresponding ttkTriangulationArray to the
   * vtkDataSet's field data. If the vtkDataSet has such a field data array,
   * then the method performs a fast check to verify that the referenced
   * triangulation is still valid for the given vtkDataSet and updates the
   * triangulation if necessary. Finally, the method returns the triangulation
   * that is referenced by the ttkTriangulationArray.
   *
   * To pass the triangulation along the pipeline, filters have to perform a
   * shallow or deep copy of an input that already has a triangulation.
   */
  ttk::Triangulation *GetTriangulation(vtkDataSet *dataSet);

  /**
   * This key can be used during the FillOutputPortInformation() call to
   * specify that an output port should produce the same data type as a
   * certain input port.
   */
  static vtkInformationIntegerKey *SAME_DATA_TYPE_AS_INPUT_PORT();

  /**
   * This method processes a pipeline request such as
   * vtkDemandDrivenPipeline::REQUEST_DATA or
   * vtkDemandDrivenPipeline::REQUEST_INFORMATION.
   *
   * It is not recommended to override this method in order to be conform
   * to the VTK/TTK pipeline model.
   */
  int ProcessRequest(vtkInformation *request,
                     vtkInformationVector **inputVectors,
                     vtkInformationVector *outputVector) override;

  /**
   * Get the output data object for a port on this algorithm.
   */
  vtkDataSet *GetOutput();
  vtkDataSet *GetOutput(int);

  /**
   * Assign a data object as input. Note that this method does not
   * establish a pipeline connection. Use SetInputConnection() to
   * setup a pipeline connection.
   */
  void SetInputData(vtkDataSet *);
  void SetInputData(int, vtkDataSet *);

  /**
   * Assign a data object as input. Note that this method does not
   * establish a pipeline connection. Use AddInputConnection() to
   * setup a pipeline connection.
   */
  void AddInputData(vtkDataSet *);
  void AddInputData(int, vtkDataSet *);

  /**
   * @brief This method tests whether the input is a nullptr.
   * If the computation is being done on multiple processes, it is possible
   * that the domain of one process or more is empty, but not others, therefore
   * in that particular case the rest of the filter will not be computed but
   * an error message will not be sent.
   *
   * @tparam inputType
   * @param input  the input to assess
   * @return int 0: error, 1: stop without error, 2: continue
   */
  template <typename inputType>
  inline int checkEmptyMPIInput(inputType *input) {
    if(!input) {
#ifdef TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI()) {
        return 1;
      } else {
#endif
        return 0;
#ifdef TTK_ENABLE_MPI
      }
#endif
    }
    return 2;
  };

protected:
  ttkAlgorithm();
  ~ttkAlgorithm() override;
  float CompactTriangulationCacheSize{0.2f};

#ifdef TTK_ENABLE_MPI
  /**
   * @brief Creates a new communicator if one of the processes doesn't contain
   * any point or cells. In this case, the RankArray is update for vertices and
   * cells to match the new ranks.
   *
   * @param input input data set
   * @return int 0 if input contains no points or no cells
   */
  int updateMPICommunicator(vtkDataSet *input);
#endif
  /**
   * This method is called in GetTriangulation, after the triangulation as been
   * created. It verifies that ghost cells and points are present and if they
   * are not, computes them.
   */

  void MPIGhostPipelinePreconditioning(vtkDataSet *input);

  /**
   * This method is called in GetTriangulation, after the triangulation as been
   * created. It verifies that several attributes necessary for MPI computation
   * are present in the pipeline and if not, computes them.
   */
  void MPIPipelinePreconditioning(vtkDataSet *input,
                                  std::vector<int> &neighbors,
                                  ttk::Triangulation *triangulation = nullptr);

  /**
   * This method checks the validity of the global identifiers given in
   * argument. A set of global identifiers is valid if the highest global id is
   * equal to the global number of simplices (without ghosts) - 1 and if the
   * lowest global id is equal to 0. Simplices can either be vertices or
   * simplices of highest dimension.
   */

  bool checkGlobalIdValidity(ttk::LongSimplexId *globalIds,
                             ttk::SimplexId simplexNumber,
                             unsigned char *ghost,
                             int *rankArray);
  /**
   * This methods generates global ids and is called during the MPI
   * preconditioning. It behaves differently for PolyData and ImageData
   * datasets.
   */

  int GenerateGlobalIds(
    vtkDataSet *input,
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> &vertGtoL,
    std::vector<int> &neighborRanks);

  /**
   * This method is called in GetTriangulation, after the triangulation as been
   * created. It retrieves several attributes from the pipeline to precondition
   * the triangulation for MPI computation.
   */

  void MPITriangulationPreconditioning(ttk::Triangulation *triangulation,
                                       vtkDataSet *input);

  /**
   * This method is called during the first pipeline pass in
   * ProcessRequest() to create empty output data objects. The data type of
   * the generated outputs is specified in FillOutputPortInformation().
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestDataObject(vtkInformation *request,
                                vtkInformationVector **inputVectors,
                                vtkInformationVector *outputVector);

  /**
   * This method is called during the second pipeline pass in
   * ProcessRequest() to provide lightweight information about the outputs
   * without any lengthy computations. For example, the data extend or the
   * number of available time steps.
   *
   * In general, it should only be necessary to override this method to
   * provide information about new vtkImageData output objects, such as
   * their extend, spacing, and origin.
   */
  virtual int
    RequestInformation(vtkInformation *ttkNotUsed(request),
                       vtkInformationVector **ttkNotUsed(inputVectors),
                       vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method is called during the third pipeline pass in
   * ProcessRequest() to update time.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int
    RequestUpdateTime(vtkInformation *ttkNotUsed(request),
                      vtkInformationVector **ttkNotUsed(inputVectors),
                      vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method is called during the fourth pipeline pass in
   * ProcessRequest() to set time dependent information.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestUpdateTimeDependentInformation(
    vtkInformation *ttkNotUsed(request),
    vtkInformationVector **ttkNotUsed(inputVectors),
    vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method is called during the fifth pipeline pass in
   * ProcessRequest() to specify which portion of its input is needed to
   * create the portion of its output that a downstream filter requested.
   *
   * In general it should not be necessary to override this method unless
   * the filter supports spatial or temporal streaming.
   */
  virtual int
    RequestUpdateExtent(vtkInformation *ttkNotUsed(request),
                        vtkInformationVector **ttkNotUsed(inputVectors),
                        vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method is called during the sixth pipeline pass in
   * ProcessRequest() to specify which outputs will currently not be
   * generated during a RequestData() call.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int
    RequestDataNotGenerated(vtkInformation *ttkNotUsed(request),
                            vtkInformationVector **ttkNotUsed(inputVectors),
                            vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method is called during the seventh pipeline pass in
   * ProcessRequest() to execute an algorithm and update the so far empty
   * output data objects.
   *
   * This method has to be overridden in order to implement the purpose of
   * the filter.
   */
  virtual int RequestData(vtkInformation *ttkNotUsed(request),
                          vtkInformationVector **ttkNotUsed(inputVectors),
                          vtkInformationVector *ttkNotUsed(outputVector)) {
    return 1;
  }

  /**
   * This method specifies the required input object data types of the
   * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
   * the port information.
   *
   * This method has to be overridden to specify the required input data
   * types.
   */
  int FillInputPortInformation(int ttkNotUsed(port),
                               vtkInformation *ttkNotUsed(info)) override {
    return 0;
  }

  /**
   * This method specifies in the port information the data type of the
   * output objects. It is possible to either explicitly specify a type by
   * adding a vtkDataObject::DATA_TYPE_NAME() key, or to pass a type of an
   * input port to an output port by adding the
   * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key.
   *
   * This method has to be overridden to specify the data types of the
   * outputs.
   */
  int FillOutputPortInformation(int ttkNotUsed(port),
                                vtkInformation *ttkNotUsed(info)) override {
    return 0;
  }
};

/// @}
