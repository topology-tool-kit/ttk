/**
 * \ingroup vtk
 * \class ttkContourAroundPoint
 * \author Christopher Kappe <kappe@cs.uni-kl.de>
 * \date January 25, 2019
 *
 * \brief TTK VTK-filter that wraps the contourAroundPoint processing package.
 *
 * VTK wrapping code for the ttk::ContourAroundPoint package.
 *
 * \param Input Input scalar field (vtkDataSet)
 * \param Output Output scalar field (vtkDataSet)
 *
 * See the related ParaView example state files for usage examples within a VTK
 * pipeline.
 *
 * \sa ttk::ContourAroundPoint
 */
#pragma once

#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ContourAroundPoint.hpp>
#include <ttkAlgorithm.h>
#include <ttkContourAroundPointModule.h> // for TTKCONTOURAROUNDPOINT_EXPORT

#include <Triangulation.h> // for tk::Triangulation::Type

class vtkInformation;
class vtkInformationVector;
class vtkUnstructuredGrid;

class TTKCONTOURAROUNDPOINT_EXPORT ttkContourAroundPoint
  : public ttkAlgorithm,
    protected ttk::ContourAroundPoint {

public:
  static ttkContourAroundPoint *New();
  vtkTypeMacro(ttkContourAroundPoint, ttkAlgorithm);

  // setter+getter macros for each parameter from the XML file
  vtkSetMacro(ui_extension, double) vtkGetMacro(ui_extension, double);
  vtkSetMacro(ui_sizeFilter, double) vtkGetMacro(ui_sizeFilter, double);
  vtkSetMacro(ui_spherical, bool) vtkGetMacro(ui_spherical, bool);

  // for the standalone (maybe unify with the above sometime)
  void SetRegionExtension(double val) {
    ui_extension = val;
  }
  void SetSizeFilter(double val) {
    ui_sizeFilter = val;
  }
  void SetSpherical(bool val) {
    ui_spherical = val;
  }

protected:
  ttkContourAroundPoint() {
    SetNumberOfInputPorts(3);
    SetNumberOfOutputPorts(2);
  }

  ~ttkContourAroundPoint() override {
  }

  // Make sure this is consistent with the XML file and the
  // `SetNumberOfInputPorts` and `SetNumberOfOutputPorts` argument
  // (used in the constructor). The return value is interpreted as
  // "PORT_REQUIREMENTS_FILLED" in vtkAlgorithm.cxx
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **iVec,
                  vtkInformationVector *oVec) override;

  /// @return Went well?
  bool preprocessFld(vtkDataSet *dataset);

  /// @return Went well?
  bool preprocessPts(vtkUnstructuredGrid *nodes, vtkUnstructuredGrid *arcs);

  /// @return Went well?
  bool process();

  /// Assemble the output object from the results of the TTK module.
  bool postprocess();

  template <typename T>
  T *getBuffer(vtkFieldData *data,
               const std::string &varName,
               int typeCode,
               const std::string &typeName) {
    auto vtkArr = data->GetAbstractArray(varName.c_str());
    const std::string dataKind
      = dynamic_cast<vtkPointData *>(data) ? "point" : "cell";
    if(!vtkArr) {
      vtkErrorMacro("The " + dataKind + "s must have data named " + varName);
      return nullptr;
    }
    if(vtkArr->GetDataType() != typeCode) {
      vtkErrorMacro(<< "The " + dataKind + " data " + varName
                         + " must be of type "
                    << typeName << " but it is "
                    << vtkArr->GetDataTypeAsString());
      return nullptr;
    }
    return reinterpret_cast<T *>(vtkArr->GetVoidPointer(0));
  }

private:
  // output region size, in percent of the maximum overlap-free size:
  // 0 --> single point, 100 --> touching neighbor regions
  double ui_extension;
  // minimum required output region size,
  // in centi-percent of the number of input field vertices:
  // 0 --> all pass, 10000 none pass
  double ui_sizeFilter;
  // name of the scalar variable of the input field
  bool ui_spherical;

  ttk::Triangulation::Type _triangTypeCode; // triangulation->getType()
  int _scalarTypeCode; // VTK type of the scalars defined on the input field
  const char *_scalarsName = nullptr;

  // referring to the input points
  std::vector<float> _coords;
  std::vector<float> _scalars;
  std::vector<float> _isovals;
  std::vector<int> _flags;

  vtkSmartPointer<vtkUnstructuredGrid> _outFld;
  vtkSmartPointer<vtkUnstructuredGrid> _outPts;

  std::vector<float> coordsBuf_{};
  std::vector<ttk::LongSimplexId> cinfosBuf_{};
  std::vector<float> scalarsBuf_{};
  std::vector<int> flagsBuf_{};
};
