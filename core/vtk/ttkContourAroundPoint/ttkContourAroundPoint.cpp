#include <ttkContourAroundPoint.h>

#include <vtkProbeFilter.h>

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

vtkStandardNewMacro(ttkContourAroundPoint)

  int ttkContourAroundPoint::doIt(std::vector<vtkDataSet *> &inputs,
                                  std::vector<vtkDataSet *> &outputs) {
  ttk::Memory memUseObj;
  _out = static_cast<vtkUnstructuredGrid *>(outputs[0]);

#if 0
  makeDummyOutput();
#else
  if(!preprocessDomain(inputs[0]))
    return 0;
  if(!preprocessConstraints(static_cast<vtkUnstructuredGrid *>(inputs[1]),
                            static_cast<vtkUnstructuredGrid *>(inputs[2])))
    return 0;
  if(!process())
    return 0;
  if(!postprocess())
    return 0;
#endif

  std::ostringstream memUseStream;
  memUseStream << std::fixed << std::setprecision(3)
               << memUseObj.getElapsedUsage() << " MB";
  dMsg(cout, memUseStream.str(), memoryMsg);
  return 1;
}

//------------------------------------------------------------------------------------------------//

bool ttkContourAroundPoint::preprocessDomain(vtkDataSet *dataset) {
  if(ui_scalars == "") {
    vtkErrorMacro(
      "A scalar variable needs to be defined on the Domain") return false;
  }

  ttk::Triangulation *triangulation
    = ttkTriangulation::getTriangulation(dataset);
  if(!triangulation) {
    vtkErrorMacro("No ttk::Triangulation could be gained from the input "
                  "field") return false;
  }

  triangulation->setWrapper(this);
  auto scalars = dataset->GetPointData()->GetAbstractArray(ui_scalars.c_str());
  assert(scalars);
  const auto errorCode
    = _wrappedModule.setupDomain(triangulation, scalars->GetVoidPointer(0));
  if(errorCode < 0) {
    vtkErrorMacro("_wrappedModule.setupDomain failed with code " << errorCode);
    return false;
  }

  _scalarTypeCode = scalars->GetDataType();
#ifndef NDEBUG
  std::ostringstream stream;
  stream << "Scalar type: " << scalars->GetDataTypeAsString() << " (code "
         << _scalarTypeCode << ")";
  dMsg(std::cout, stream.str().c_str());
#endif

  double *bounds = dataset->GetBounds();
  double size = 1.;
  for(std::size_t dim = 0; dim < 3; ++dim) {
    const double range = bounds[dim * 2 + 1] - bounds[dim * 2];
    size *= range == 0 ? 1 : range;
  }
  _domainBbSize = size;

  return true;
}

//------------------------------------------------------------------------------------------------//

bool ttkContourAroundPoint::preprocessConstraints(vtkUnstructuredGrid *nodes,
                                                  vtkUnstructuredGrid *arcs) {
  // ---- Point data ---- //

  auto points = nodes->GetPoints();
  if(points->GetDataType() != VTK_FLOAT) {
    vtkErrorMacro("The point coordinates must be of type float") return false;
  }
  auto coords = reinterpret_cast<float *>(points->GetData()->GetVoidPointer(0));

  auto pData = nodes->GetPointData();
  //  auto sizeBuf = getBuffer<float>(pData, "RegionSize", VTK_FLOAT, "float");
  auto sizeBuf = getBuffer<int>(
    pData, "RegionSize", VTK_INT, "int"); // TODO Document precise meaning
  auto codeBuf = getBuffer<int>(pData, "CriticalType", VTK_INT, "int");
  auto scalarBuf = getBuffer<float>(pData, "Scalar", VTK_FLOAT, "float");
  if(!sizeBuf || !codeBuf || !scalarBuf)
    return false;

  // ---- Cell data ---- //

  auto cells = arcs->GetCells();
  const auto maxNvPerC
    = cells->GetMaxCellSize(); // NOTE There is no equivalent for min :-(
  if(maxNvPerC != 2) {
    vtkErrorMacro(
      "The points must come in pairs but there is at least one cell with "
      + std::to_string(maxNvPerC) + " points") return false;
  }

  auto cData = arcs->GetCellData();
  auto c2p = getBuffer<int>(cData, "upNodeId", VTK_INT, "int");
  auto c2q = getBuffer<int>(cData, "downNodeId", VTK_INT, "int");
  if(!c2p || !c2q)
    return false;

  // ---- Loop over pairs ---- //

  static constexpr int minCode = 0;
  static constexpr int maxCode = 3;
  const double sadFac
    = ui_extension * 0.01; // factor for the saddle or mean of min and max
  const double extFac = 1 - sadFac; // factor for the extreme point
  const double minSize = _domainBbSize * ui_sizeFilter * 0.0001;

  _isovals.resize(0);
  _coords.resize(0);
  _flags.resize(0);

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
      if(sizeBuf[ext] < minSize)
        continue;
      _isovals.push_back(scalarBuf[ext] * extFac + scalarBuf[sad] * sadFac);
      const auto point = &coords[ext * 3];
      _coords.push_back(point[0]);
      _coords.push_back(point[1]);
      _coords.push_back(point[2]);
      _flags.push_back(codeBuf[ext] == minCode ? 0 : 1);
    } else // min-max pair
    {
      const auto pVal = scalarBuf[p];
      const auto qVal = scalarBuf[q];
      vtkWarningMacro(<< "Arc " << c << " joins a minimum and a maximum")
        const auto cVal
        = (pVal + qVal) / 2;
      if(sizeBuf[p] >= minSize) {
        _isovals.push_back(pVal * extFac + cVal * sadFac);
        const auto point = &coords[p * 3];
        _coords.push_back(point[0]);
        _coords.push_back(point[1]);
        _coords.push_back(point[2]);
        _flags.push_back(pCode == minCode ? 0 : 1);
      }
      if(sizeBuf[q] >= minSize) {
        _isovals.push_back(qVal * extFac + cVal * sadFac);
        const auto point = &coords[q * 3];
        _coords.push_back(point[0]);
        _coords.push_back(point[1]);
        _coords.push_back(point[2]);
        _flags.push_back(qCode == minCode ? 0 : 1);
      }
    }
  }

  const auto errorCode = _wrappedModule.setupConstraints(
    _coords.data(), _isovals.data(), _isovals.size(), _flags.data());
  if(errorCode < 0) {
    vtkErrorMacro("setInputPoints failed with code " << errorCode);
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------------------------//

bool ttkContourAroundPoint::process() {
  _wrappedModule.setWrapper(this);
  int errorCode = 0; // In TTK, negative is bad.
  switch(_scalarTypeCode) {
    ttkTemplateMacro((errorCode = _wrappedModule.execute<VTK_TT>()));
  }
  if(errorCode < 0) {
    vtkErrorMacro("_wrappedModule.execute failed with code "
                  << errorCode) return false;
  }

  return true;
}

//------------------------------------------------------------------------------------------------//

bool ttkContourAroundPoint::postprocess() {
  float *coordsBuf;
  ttk::SimplexId nv;
  ttk::SimplexId *cinfosBuf;
  ttk::SimplexId nc;
  float *scalarsBuf;
  int *flagsBuf;
  _wrappedModule.getOutputField(
    coordsBuf, nv, cinfosBuf, nc, scalarsBuf, flagsBuf);
  if(nv == 0) // very fine area filter
    return true;

  const auto nvPerC = cinfosBuf[0];
  if(!(nvPerC == 2 || nvPerC == 3)) {
    vtkErrorMacro("Invalid number of vertices per cell: "
                  + std::to_string(nvPerC) + ". Must be 2 or 3") return false;
  }
  const auto cellType = nvPerC == 2 ? VTK_LINE : VTK_TRIANGLE;

#ifndef NDEBUG
  if(vtkSmartPointer<vtkPoints>::New()->GetDataType() != VTK_FLOAT) {
    vtkErrorMacro("The API has changed! We have expected the default "
                  "coordinate type to be float");
    return false;
  }

  for(ttk::SimplexId c = 0, i = 0; c < nc; ++c) {
    const auto n = cinfosBuf[i];
    if(n != nvPerC) {
      vtkErrorMacro("Output cell " + std::to_string(c) + " has "
                    + std::to_string(n) + " vertices. " + "Expected "
                    + std::to_string(nvPerC));
      return false;
    }
    i += n + 1;
  }
#endif

  const std::size_t nEntriesPerC = nvPerC + 1;
  auto cinfosBufVtk = reinterpret_cast<vtkIdType *>(cinfosBuf);
  if(!std::is_same<ttk::SimplexId, vtkIdType>::value) {
    const std::size_t n = nc * nEntriesPerC;
    cinfosBufVtk = new vtkIdType[n];
    for(std::size_t i = 0; i < n; ++i)
      cinfosBufVtk[i] = vtkIdType(cinfosBuf[i]);
    delete[] cinfosBuf;
  }

  // Pass ownership of the heap-allocated raw array to the respective
  // vtkDataArray.
  const int wantSave = 0;
  // Use `delete[]` instead of the VTK default `free()`.
  // (The enum is independent of the template type - just use float.)
  const int delMethod
    = VTK_AOS_DATA_ARRAY_TEMPLATE<float>::DeleteMethod::VTK_DATA_ARRAY_DELETE;

  auto points = vtkSmartPointer<vtkPoints>::New();
  auto coordArr = vtkSmartPointer<vtkFloatArray>::New();
  coordArr->SetNumberOfComponents(3);
  coordArr->SetArray(coordsBuf, nv * 3, wantSave, delMethod);
  points->SetData(coordArr);
  _out->SetPoints(points);

  auto cells = vtkSmartPointer<vtkCellArray>::New();
  auto cinfoArr = vtkSmartPointer<vtkIdTypeArray>::New();
  cinfoArr->SetArray(cinfosBufVtk, nc * nEntriesPerC, wantSave, delMethod);
  cells->SetCells(nc, cinfoArr);
  _out->SetCells(cellType, cells);

  auto scalarArr = vtkFloatArray::New();
  scalarArr->SetArray(scalarsBuf, nv, wantSave, delMethod);
  scalarArr->SetName(ui_scalars.c_str());
  _out->GetPointData()->AddArray(scalarArr);

  auto flagArr = vtkIntArray::New();
  flagArr->SetArray(flagsBuf, nv, wantSave, delMethod);
  flagArr->SetName("isMax");
  _out->GetPointData()->AddArray(flagArr);

  return true;
}

//------------------------------------------------------------------------------------------------//

void ttkContourAroundPoint::makeDummyOutput() {
  auto points = vtkSmartPointer<vtkPoints>::New();
  // longitude (x) in [0,360), latitude (y) in (-90,90) (the poles are
  // singularites)
  points->InsertNextPoint(180, 45, 0); // north center
  points->InsertNextPoint(90, -45, 0); // south east
  points->InsertNextPoint(270, -45, 0); // south west
  _out->SetPoints(points);

  auto cells = vtkSmartPointer<vtkCellArray>::New();
  // Anonymous arrays (passed directly to the function) would only be possible
  // if `InsertNextCell` were overloaded with an
  // std::initializer_list<vtkIdType> argument :-(
  const vtkIdType c0[] = {0, 1};
  cells->InsertNextCell(2, c0);
  const vtkIdType c1[] = {1, 2};
  cells->InsertNextCell(2, c1);
  const vtkIdType c2[] = {2, 0};
  cells->InsertNextCell(2, c2);
  _out->SetCells(VTK_LINE, cells);

  auto scalarArr = vtkFloatArray::New();
  static constexpr float placeholder = 0.1337;
  for(int i = 0; i < 3; ++i)
    scalarArr->InsertNextValue(placeholder);
  scalarArr->SetName(ui_scalars.c_str());
  _out->GetPointData()->AddArray(scalarArr);
}
