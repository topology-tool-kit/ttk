/// \ingroup vtk
/// \class ttkAlgorithm
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2019.
///
/// \brief Baseclass of all VTK filters that wrap ttk modules.
///
/// This is an abstract vtkAlgorithm that provides standardized input/output
/// managment for VTK wrappers of ttk filters. The class also provides a static
/// method to retrieve a ttk::Triangulation of a vtkDataSet.

#pragma once

// VTK Module
#include <ttkAlgorithmModule.h>

// std includes
#include <unordered_map>

// VTK Includes
#include <vtkAlgorithm.h>
class vtkCellArray;
class vtkCommand;
class vtkDataSet;
class vtkInformation;
class vtkInformationIntegerKey;
class vtkPoints;

template <class d0>
class vtkSmartPointer;

// Base Includes
#include <Debug.h>

namespace ttk {
  class Triangulation;
}
// #include <Triangulation.h>

class TTKALGORITHM_EXPORT ttkAlgorithm : public vtkAlgorithm,
                                         virtual public ttk::Debug {
private:
  /**
   * A static registry that maps owners (e.g. vtkCellArrays or vtkImageData
   * objects) to a ttk::Triangulation object. The registry also checks if a
   * triangulation needs to be updated in case the owner was modified since
   * initialization, and it also stores an event listener that automatically
   * deletes a triangulation if its corresponding owner is deleted.
   */
  static std::unordered_map<void *,
                            std::tuple<ttk::Triangulation,
                                       vtkObject *,
                                       vtkSmartPointer<vtkCommand>,
                                       vtkMTimeType>>
    DataSetToTriangulationMap;

  int ThreadNumber{1};
  bool UseAllCores{true};

  /**
   * This function checks if the registry contains a triangulation for a given
   * owner. It also checks if the triangulation would need an update, in which
   * case the triangulation and its auxiliary objects are deleted from the
   * registry (now a new triangulation can be recreated from scratch).
   */
  ttk::Triangulation *FindTriangulation(void *key);

  /**
   * This function creates a ttk::Triangulation object for a given.
   * Specifically, it initializes either an explicit ttk::Triangulation (in case
   * points and cells are provided), or an implicit triangulation (in case owner
   * is a vtkImageData object).
   */
  ttk::Triangulation *InitTriangulation(void *key,
                                        vtkObject *owner,
                                        vtkPoints *points = nullptr,
                                        vtkCellArray *cells = nullptr);

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
   * (overriden by UseAllCores member).
   */
  void SetThreadNumber(int threadNumber) {
    this->ThreadNumber = threadNumber;
    this->UpdateThreadNumber();
  };

  /**
   * Controls if the base code should use all available cores.
   */
  void SetUseAllCores(bool useAllCores) {
    this->UseAllCores = useAllCores;
    this->UpdateThreadNumber();
  };

  /**
   * Controls the debug level used by algorithms that are invoked by the VTK
   * wrapper.
   */
  void SetDebugLevel(int debugLevel) {
    this->setDebugLevel(debugLevel); // from ttk::Debug
    this->Modified();
  };

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
                                 vtkInformationVector **inputVectors,
                                 const int &inputPort = 0);

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
  ttk::Triangulation *GetTriangulation(vtkDataSet *object);

  /**
   * This key can be used during the FillOutputPortInfomration() call to
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

protected:
  ttkAlgorithm();
  virtual ~ttkAlgorithm();

  /**
   * This method is called during the first pipeline pass in
   * ProcessRequest() to create empty output data objects. The data type of
   * the generated outputs is specified in FillOutputPortInfomration().
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
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVectors,
                                 vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the third pipeline pass in
   * ProcessRequest() to update time.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestUpdateTime(vtkInformation *request,
                                vtkInformationVector **inputVectors,
                                vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the fourth pipeline pass in
   * ProcessRequest() to set time dependent information.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int
    RequestUpdateTimeDependentInformation(vtkInformation *request,
                                          vtkInformationVector **inputVectors,
                                          vtkInformationVector *outputVector) {
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
  virtual int RequestUpdateExtent(vtkInformation *request,
                                  vtkInformationVector **inputVectors,
                                  vtkInformationVector *outputVector) {
    return 1;
  };

  /**
   * This method is called during the sixth pipeline pass in
   * ProcessRequest() to specify which outputs will currently not be
   * generated during a RequestData() call.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestDataNotGenerated(vtkInformation *request,
                                      vtkInformationVector **inputVectors,
                                      vtkInformationVector *outputVector) {
    return 1;
  };

  /**
   * This method is called during the seventh pipeline pass in
   * ProcessRequest() to execute an algorithm and update the so far empty
   * output data objects.
   *
   * This method has to be overridden in order to implement the purpose of
   * the filter.
   */
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVectors,
                          vtkInformationVector *outputVector) {
    return 1;
  };

  /**
   * This method specifies the required input object data types of the
   * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
   * the port information.
   *
   * This method has to be overridden to specify the required input data
   * types.
   */
  virtual int FillInputPortInformation(int port,
                                       vtkInformation *info) override {
    return 0;
  };

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
  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override {
    return 0;
  };
};
