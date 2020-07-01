#include <ttkContourAroundPoint.h>

#include <vtkVersion.h>

#include <vtkProbeFilter.h> // TODO can be removed?

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

#if VTK_MAJOR_VERSION >= 7
#include <vtkAOSDataArrayTemplate.h>
#define VTK_AOS_DATA_ARRAY_TEMPLATE vtkAOSDataArrayTemplate
#else
#include <vtkDataArrayTemplate.h>
#define VTK_AOS_DATA_ARRAY_TEMPLATE vtkDataArrayTemplate
#endif

#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>

#include <vtkSmartPointer.h>

#include <cassert>
#include <type_traits>

vtkStandardNewMacro(ttkContourAroundPoint);

int ttkContourAroundPoint::doIt(std::vector<vtkDataSet *> &inputs,
                                std::vector<vtkDataSet *> &outputs) {
  ttk::Memory memUseObj;
  _outFld = static_cast<vtkUnstructuredGrid *>(outputs[0]);
  _outPts = static_cast<vtkUnstructuredGrid *>(outputs[1]);

  if(!preprocessFld(inputs[0]))
    return 0;
  if(!preprocessPts(static_cast<vtkUnstructuredGrid *>(inputs[1]),
                    static_cast<vtkUnstructuredGrid *>(inputs[2])))
    return 0;
  if(!process())
    return 0;
  if(!postprocess())
    return 0;

  std::ostringstream memUseStream;
  memUseStream << std::fixed << std::setprecision(3)
               << memUseObj.getElapsedUsage() << " MB";
  dMsg(cout, memUseStream.str(), memoryMsg);
  return 1;
}

//----------------------------------------------------------------------------//

bool ttkContourAroundPoint::preprocessFld(vtkDataSet *dataset) {
  if(ui_scalars == "") {
    vtkErrorMacro("A scalar variable needs to be defined on the Domain");
    return false;
  }

  ttk::Triangulation *triangulation
    = ttkTriangulation::getTriangulation(dataset);
  if(!triangulation) {
    vtkErrorMacro("No ttk::Triangulation could be gained from the input field");
    return false;
  }

  triangulation->setWrapper(this);
  auto scalars = dataset->GetPointData()->GetAbstractArray(ui_scalars.c_str());
  assert(scalars);
  const double radius = ui_spherical ? -1. : 0.;
  const auto errorCode = _wrappedModule.setInputField(
        triangulation, scalars->GetVoidPointer(0), ui_sizeFilter, radius);
  if(errorCode < 0) {
    vtkErrorMacro("_wrappedModule.setInputField failed with code " << errorCode);
    return false;
  }

  _scalarTypeCode = scalars->GetDataType();
  std::ostringstream stream;
  stream << "Scalar type: " << scalars->GetDataTypeAsString() << " (code "
         << _scalarTypeCode << ")";
  dMsg(std::cout, stream.str().c_str(), detailedInfoMsg);
  return true;
}

//----------------------------------------------------------------------------//

bool ttkContourAroundPoint::preprocessPts(vtkUnstructuredGrid *nodes,
                                          vtkUnstructuredGrid *arcs) {
  // ---- Point data ---- //

  auto points = nodes->GetPoints();
  if(points->GetDataType() != VTK_FLOAT) {
    vtkErrorMacro("The point coordinates must be of type float");
    return false;
  }
  auto coords = reinterpret_cast<float *>(points->GetData()->GetVoidPointer(0));

  auto pData = nodes->GetPointData();
  auto scalarBuf = getBuffer<float>(pData, "Scalar", VTK_FLOAT, "float");
  auto codeBuf = getBuffer<int>(pData, "CriticalType", VTK_INT, "int");
  if(!scalarBuf || !codeBuf)
    return false;

  // ---- Cell data ---- //

#ifndef NDEBUG // each arc should of course be defined by exactly two vertices
  auto cells = arcs->GetCells();
  const auto maxNvPerC = cells->GetMaxCellSize();
  if(maxNvPerC != 2) {
    vtkErrorMacro(
      "The points must come in pairs but there is at least one cell with "
      + std::to_string(maxNvPerC) + " points");
    return false;
  }
  // TODO Check for minNvPerC != 2
#endif

  auto cData = arcs->GetCellData();
  auto c2p = getBuffer<int>(cData, "upNodeId", VTK_INT, "int");
  auto c2q = getBuffer<int>(cData, "downNodeId", VTK_INT, "int");
  if(!c2p || !c2q)
    return false;

  // ---- Loop over pairs ---- //

  static constexpr int minCode = 0;
  static constexpr int maxCode = 3;
  const double sadFac = ui_extension * 0.01; // saddle or mean of min and max
  const double extFac = 1 - sadFac; // factor for the extreme point

  _coords.resize(0);
  _scalars.resize(0);
  _isovals.resize(0);
  _flags.resize(0);
  
  auto addPoint =
      [this, coords, scalarBuf, codeBuf](int p, float isoval, int code) {
    const auto point = &coords[p * 3];
    _coords.push_back(point[0]);
    _coords.push_back(point[1]);
    _coords.push_back(point[2]);
    _scalars.push_back(scalarBuf[p]);
    _isovals.push_back(isoval);
    _flags.push_back(code == minCode ? 0 : 1);
  };
  
  const vtkIdType nc = arcs->GetNumberOfCells();
  for(vtkIdType c = 0; c < nc; ++c) {
    const auto p = c2p[c];
    const auto q = c2q[c];
    const auto pCode = codeBuf[p];
    const auto qCode = codeBuf[q];

    const bool pIsSad = pCode != minCode && pCode != maxCode;
    const bool qIsSad = qCode != minCode && qCode != maxCode;
    if(pIsSad || qIsSad) {
      if(pIsSad && qIsSad) // two saddles
        continue;
      // extremum and saddle
      const auto ext = pIsSad ? q : p;
      const auto sad = pIsSad ? p : q;
      const float isoval = scalarBuf[ext] * extFac + scalarBuf[sad] * sadFac;
      addPoint(ext, isoval, codeBuf[ext]);
    } else { // min-max pair
      vtkWarningMacro(<< "Arc " << c << " joins a minimum and a maximum");
      const auto pVal = scalarBuf[p];
      const auto qVal = scalarBuf[q];
      const auto cVal = (pVal + qVal) / 2;
      addPoint(p, pVal * extFac + cVal * sadFac, pCode);
      addPoint(q, qVal * extFac + cVal * sadFac, qCode);
    }
  }

  auto np = _scalars.size();
  const auto errorCode = _wrappedModule.setInputPoints(
    _coords.data(), _scalars.data(), _isovals.data(), _flags.data(), np);
  if(errorCode < 0) {
    vtkErrorMacro("setInputPoints failed with code " << errorCode);
    return false;
  }
  return true;
}

//----------------------------------------------------------------------------//

bool ttkContourAroundPoint::process() {
  _wrappedModule.setWrapper(this);
  int errorCode = 0; // In TTK, negative is bad.
  switch(_scalarTypeCode) {
    vtkTemplateMacro((errorCode = _wrappedModule.execute<VTK_TT>()));
  }
  if(errorCode < 0) {
    vtkErrorMacro("_wrappedModule.execute failed with code " << errorCode);
    return false;
  }
  return true;
}

//----------------------------------------------------------------------------//

bool ttkContourAroundPoint::postprocess() {
  ttk::SimplexId *cinfosBuf;
  ttk::SimplexId nc;
  float *coordsBuf;
  float *scalarsBuf;
  int *flagsBuf;
  ttk::SimplexId nv;
  _wrappedModule.getOutputContours(
    cinfosBuf, nc, coordsBuf, scalarsBuf, flagsBuf, nv);
  if(nc == 0) // very fine area filter
    return true;

  // Pass ownership of the heap-allocated raw array to the respective
  // vtkDataArray.
  const int wantSave = 0;
  // Use `delete[]` instead of the VTK default `free()`.
  // (The enum is independent of the template type - just use float.)
  const int delMethod
    = VTK_AOS_DATA_ARRAY_TEMPLATE<float>::DeleteMethod::VTK_DATA_ARRAY_DELETE;

  // ---- Cell data (output 0) ---- //

  int *ctypes = new int[nc];

  vtkIdType cinfoCounter = 0;
  for(std::size_t c = 0; c < nc; ++c) {
    const auto nvOfCell = cinfosBuf[cinfoCounter];
    assert(nvOfCell >= 2 && nvOfCell <= 3); // ensured in _wrappedModule
    ctypes[c] = nvOfCell == 2 ? VTK_LINE : VTK_TRIANGLE;
    cinfoCounter += nvOfCell + 1;
  }
  const vtkIdType cinfosSize = cinfoCounter;

  auto cinfosBufVtk = reinterpret_cast<vtkIdType *>(cinfosBuf);
  if(!std::is_same<ttk::SimplexId, vtkIdType>::value) { // unlikely
    // Actually a warning would be in order-
    // what if conversion is not possible (e.g. too large indices)?
    cinfosBufVtk = new vtkIdType[cinfosSize];
    for(std::size_t i = 0; i < cinfosSize; ++i)
      cinfosBufVtk[i] = vtkIdType(cinfosBuf[i]);
    delete[] cinfosBuf;
  }

  auto cells = vtkSmartPointer<vtkCellArray>::New();
  auto cinfoArr = vtkSmartPointer<vtkIdTypeArray>::New();
  cinfoArr->SetArray(cinfosBufVtk, cinfosSize, wantSave, delMethod);
  cells->SetCells(nc, cinfoArr);
  _outFld->SetCells(ctypes, cells);

  // ---- Point data (output 0) ---- //

  if(vtkSmartPointer<vtkPoints>::New()->GetDataType() != VTK_FLOAT) {
    vtkErrorMacro("The API has changed! We have expected the default "
                  "coordinate type to be float");
    return false;
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  auto coordArr = vtkSmartPointer<vtkFloatArray>::New();
  coordArr->SetNumberOfComponents(3);
  coordArr->SetArray(coordsBuf, nv * 3, wantSave, delMethod);
  points->SetData(coordArr);
  _outFld->SetPoints(points);

  auto scalarArr = vtkFloatArray::New();
  scalarArr->SetArray(scalarsBuf, nv, wantSave, delMethod);
  scalarArr->SetName(ui_scalars.c_str());
  _outFld->GetPointData()->AddArray(scalarArr);

  auto flagArr = vtkIntArray::New();
  flagArr->SetArray(flagsBuf, nv, wantSave, delMethod);
  flagArr->SetName("isMax");
  _outFld->GetPointData()->AddArray(flagArr);

  // ---- Output 1 (added in a later revision of the algo) ---- //

  // re-using the variables from above
  _wrappedModule.getOutputCentroids(coordsBuf, scalarsBuf, flagsBuf, nv);

  points = vtkSmartPointer<vtkPoints>::New();
  coordArr = vtkSmartPointer<vtkFloatArray>::New();
  coordArr->SetNumberOfComponents(3);
  coordArr->SetArray(coordsBuf, nv * 3, wantSave, delMethod);
  points->SetData(coordArr);
  _outPts->SetPoints(points);

  scalarArr = vtkFloatArray::New();
  scalarArr->SetArray(scalarsBuf, nv, wantSave, delMethod);
  scalarArr->SetName(ui_scalars.c_str());
  _outPts->GetPointData()->AddArray(scalarArr);

  flagArr = vtkIntArray::New();
  flagArr->SetArray(flagsBuf, nv, wantSave, delMethod);
  flagArr->SetName("isMax");
  _outPts->GetPointData()->AddArray(flagArr);

  return true;
}
