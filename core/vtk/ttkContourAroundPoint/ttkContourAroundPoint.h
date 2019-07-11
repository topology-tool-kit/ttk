/**
 * \ingroup vtk
 * \class ttkContourAroundPoint
 * \author Christopher Kappe <kappe@cs.uni-kl.de>
 * \date January 25, 2019
 *
 * \brief TTK VTK-filter that wraps the contourAroundPoint processing package.
 *
 * VTK wrapping code for the @ContourAroundPoint package.
 *
 * \param Input Input scalar field (vtkDataSet)
 * \param Output Output scalar field (vtkDataSet)
 *
 * This filter can be used as any other VTK filter (for instance, by using the
 * sequence of calls SetInputData(), Update(), GetOutput()).
 *
 * See the related ParaView example state files for usage examples within a VTK
 * pipeline.
 *
 * \sa ttk::ContourAroundPoint
 */
#pragma once

// VTK includes
#include <vtkDataSetAlgorithm.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

// TTK includes
#include <ContourAroundPoint.hpp>
#include <ttkWrapper.h>

// See the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkContourAroundPoint
#else
class ttkContourAroundPoint
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkContourAroundPoint *New();
  vtkTypeMacro(ttkContourAroundPoint, vtkDataSetAlgorithm)

    // BEGIN default ttk setters
    vtkSetMacro(debugLevel_, int)

      void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // END default ttk setters

  // BEGIN set-getters macros for each parameter from the
  // ServerManagerConfiguration XML file.
  vtkSetMacro(ui_sizeFilter, double) vtkGetMacro(ui_sizeFilter, double)

    vtkSetMacro(ui_extension, double) vtkGetMacro(ui_extension, double)

      vtkSetMacro(ui_scalars, std::string) vtkGetMacro(ui_scalars, std::string)
    // END set-getters macros

    // By default, this filter has one input and one output, of the same type.
    // Here, you can override the number of input/output ports and their
    // respective type. Make sure this is consistent with the
    // ServerManagerConfiguration XML file and the
    // `SetNumberOfInputPorts`/`SetNumberOfOutputPorts` argument (used in the
    // constructor). The return value is interpreted as
    // "PORT_REQUIREMENTS_FILLED" in vtkAlgorithm.cxx
    int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        return 1;
      case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        return 1;
      case 2:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        return 1;
    }
    return 0;
  }
  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        return 1;
    }
    return 0;
  }

  /// Override this method in order to always prepend a class- and
  /// debugLevel-specific prefix and include a line end after the message.
  virtual int dMsg(std::ostream &stream,
                   std::string msg,
                   const int &debugLevel = infoMsg) const override {
    if(debugLevel > debugLevel_ && debugLevel > ttk::globalDebugLevel_)
      return 0;
    stream << "[ttkContourAroundPoint] ";
    switch(debugLevel) {
      case fatalMsg:
        stream << "Error: ";
        break; // something went wrong
      case timeMsg:
        stream << "Time consumption: ";
        break; // x.yyy s
      case memoryMsg:
        stream << "Memory usage: ";
        break; // x.yyy MB
    }
    const auto res = ttk::Debug::dMsg(stream, msg, debugLevel);
    stream << std::endl;
    return res;
  }

protected:
  ttkContourAroundPoint() {
    UseAllCores = true;
    ThreadNumber = 1;
    debugLevel_ = 3;

    SetNumberOfInputPorts(3);
    SetNumberOfOutputPorts(1);
  }

  ~ttkContourAroundPoint() {
  }

  TTK_SETUP();

  /// @return Went well?
  bool preprocessDomain(vtkDataSet *dataset);

  /// @return Went well?
  bool preprocessConstraints(vtkUnstructuredGrid *nodes,
                             vtkUnstructuredGrid *arcs);

  /// @return Went well?
  bool process();

  /// Assemble the output object from the results of the TTK module.
  bool postprocess();

  /// For testing the general pipeline stuff, without having to execute the
  /// wrapped module.
  void makeDummyOutput();

  template <typename T>
  T *getBuffer(vtkFieldData *data,
               const std::string &varName,
               int typeCode,
               const std::string &typeName) {
    auto vtkArr = data->GetAbstractArray(varName.c_str());
    const std::string dataKind
      = dynamic_cast<vtkPointData *>(data) ? "point" : "cell";
    if(!vtkArr) {
      vtkErrorMacro("The " + dataKind + "s must have data named "
                    + varName) return nullptr;
    }
    if(vtkArr->GetDataType() != typeCode) {
      vtkErrorMacro(<< "The " + dataKind + " data " + varName
                         + " must be of type "
                    << typeName << " but it is "
                    << vtkArr->GetDataTypeAsString()) return nullptr;
    }
    return reinterpret_cast<T *>(vtkArr->GetVoidPointer(0));
  }

private:
  double ui_sizeFilter; // 10000 [c%] --> keep everything
  double ui_extension; // 100 [%] --> maximum contour size
  std::string ui_scalars; // name of the scalar variable of the Domain

  int _scalarTypeCode; // VTK type of the scalars defined on the Domain
  double _domainBbSize; // size of the bounding box of the domain
  std::vector<float> _coords;
  std::vector<float> _isovals;
  std::vector<int> _flags;
  vtkSmartPointer<vtkUnstructuredGrid> _out;

  ttk::ContourAroundPoint _wrappedModule;
};
