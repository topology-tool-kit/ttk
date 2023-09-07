/// \ingroup vtk
/// \class ttkDistanceMatrixDistortion
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::DistanceMatrixDistortion module.
///
/// This VTK filter uses the ttk::DistanceMatrixDistortion module to compute the
/// distortion between two distance matrices representing the same points (for
/// instance in low and high dimensions), according to the SIM formula.
///
/// \param Input0 vtkTable.
/// \param Input1 vtkTable.
/// \param Output vtkTable.
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
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DistanceMatrixDistortion
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkDistanceMatrixDistortionModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <DistanceMatrixDistortion.h>

class TTKDISTANCEMATRIXDISTORTION_EXPORT ttkDistanceMatrixDistortion
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::DistanceMatrixDistortion // and we inherit from the base
                                            // class
{
private:
  // To store which column names we want to extract from the inputs.
  std::vector<std::string> ScalarFieldsHigh{}, ScalarFieldsLow{};
  std::string RegexpStringHigh{".*"}, RegexpStringLow{".*"};
  bool SelectFieldsWithRegexpHigh{false};
  bool SelectFieldsWithRegexpLow{false};

public:
  // Variable to choose which columns to use. Two different inputs:
  // the low one and the high one, hence duplicate variables.
  vtkSetMacro(SelectFieldsWithRegexpHigh, bool);
  vtkGetMacro(SelectFieldsWithRegexpHigh, bool);
  vtkSetMacro(SelectFieldsWithRegexpLow, bool);
  vtkGetMacro(SelectFieldsWithRegexpLow, bool);
  vtkSetMacro(DoNotNormalize, bool);
  vtkGetMacro(DoNotNormalize, bool);

  vtkSetMacro(RegexpStringHigh, const std::string &);
  vtkGetMacro(RegexpStringHigh, std::string);
  vtkSetMacro(RegexpStringLow, const std::string &);
  vtkGetMacro(RegexpStringLow, std::string);

  // Two functions because two inputs.
  void ClearScalarFieldsHigh() {
    ScalarFieldsHigh.clear();
    Modified();
  }
  void ClearScalarFieldsLow() {
    ScalarFieldsLow.clear();
    Modified();
  }

  static ttkDistanceMatrixDistortion *New();
  vtkTypeMacro(ttkDistanceMatrixDistortion, ttkAlgorithm);

  // Two functions because two inputs. Used for columns selection.
  void SetScalarFieldsHigh(const std::string &s) {
    ScalarFieldsHigh.emplace_back(s);
    Modified();
  }
  void SetScalarFieldsLow(const std::string &s) {
    ScalarFieldsLow.emplace_back(s);
    Modified();
  }

protected:
  ttkDistanceMatrixDistortion();
  ~ttkDistanceMatrixDistortion() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
