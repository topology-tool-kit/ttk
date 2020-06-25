#include <ttkTriangulation.h>
#include <ttkTriangulationAlgorithm.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

/**
 * @brief Auxiliary struct to check if VTK cells are simplices.
 */
struct CellChecker : public vtkObject
// inherit to be able to use vtkWarningMacro
{
  vtkTypeMacro(CellChecker, vtkObject);

  CellChecker() = default;
  static CellChecker *New() {
    return new CellChecker;
  }

  static std::string
    getErrMsg() // in C++ > 11: static constexpr errMsg = "foo";
  {
    return "Input explicit mesh is not (completely) simplicial. "
           "Please apply the Tetrahedralize filter before proceeding with TTK. "
           "The application behavior may be unspecified otherwise.";
  }

  int check1d(vtkDataSet *ds) {
    // e.g. VTK_POLY_LINE is not supported (as from vtkPolyData::GetLines)
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c) {
      if(ds->GetCellType(c) != VTK_LINE) {
        vtkWarningMacro(<< getErrMsg());
        return -1;
      }
    }
    return 0;
  }

  int check2d(vtkDataSet *ds) {
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c) {
      const auto type = ds->GetCellType(c);
      if(type != VTK_TRIANGLE && type != VTK_LINE) {
        vtkWarningMacro(<< getErrMsg());
        return -1;
      }
    }
    return 0;
  }

  int check3d(vtkDataSet *ds) {
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c) {
      const auto type = ds->GetCellType(c);
      if(type != VTK_TETRA && type != VTK_TRIANGLE && type != VTK_LINE) {
        vtkWarningMacro(<< getErrMsg());
        return -1;
      }
    }
    return 0;
  }
};

ttkTriangulation::ttkTriangulation() {

  inputDataSet_ = nullptr;
  hasAllocated_ = false;
  triangulation_ = nullptr;
}

ttkTriangulation::~ttkTriangulation() {

  if((triangulation_) && (hasAllocated_))
    delete triangulation_;
}

int ttkTriangulation::allocate() {

  triangulation_ = new Triangulation;
  triangulation_->setDebugLevel(0);
  hasAllocated_ = true;

  return 0;
}

int ttkTriangulation::deepCopy(vtkDataObject *other) {

  if((triangulation_) && (hasAllocated_)) {
    delete triangulation_;
  }

  allocate();

  ttkTriangulation *otherTriangulation
    = dynamic_cast<ttkTriangulation *>(other);

  if(otherTriangulation) {
    // the other object has already been enhanced
    // copy the triangulation object from the other
    (*triangulation_) = *(otherTriangulation->triangulation_);
  } else {
    triangulation_ = nullptr;
  }

  // populate the data-structure
  if((triangulation_) && (triangulation_->isEmpty())) {
    // populate the triangulation data-structure
    setInputData((vtkDataSet *)other);
  }

  return 0;
}

Triangulation *ttkTriangulation::getTriangulation(vtkDataSet *other) {

  if(!other)
    return NULL;

  ttkTriangulation *otherTriangulation
    = dynamic_cast<ttkTriangulation *>(other);

  if(otherTriangulation) {
    return otherTriangulation->getTriangulation();
  } else {
    Debug d;
    stringstream msg;
    msg << "[ttkTriangulation] Unsupported input VTK class `"
        << other->GetClassName() << "' (ref=" << other->GetDataObjectType()
        << ")" << endl;
    d.dMsg(cerr, msg.str(), Debug::fatalMsg);
  }

  return NULL;
}

vtkUnstructuredGrid *ttkTriangulation::getVtkUnstructuredGrid() {

  if(!triangulation_)
    return NULL;

  if((inputDataSet_->GetDataObjectType() == VTK_IMAGE_DATA)
     || (inputDataSet_->GetDataObjectType() == VTK_POLY_DATA)) {

    Timer t;

    vtkUnstructuredGrid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkPoints_ = vtkSmartPointer<vtkPoints>::New();

    // here we need to create an explicit representation of the triangulation
    SimplexId vertexNumber = triangulation_->getNumberOfVertices();
    vtkPoints_->SetNumberOfPoints(vertexNumber);

    float p[3];
    for(SimplexId i = 0; i < vertexNumber; i++) {
      triangulation_->getVertexPoint(i, p[0], p[1], p[2]);
      vtkPoints_->SetPoint(i, p[0], p[1], p[2]);
    }

    vtkUnstructuredGrid_->SetPoints(vtkPoints_);

    SimplexId cellNumber = triangulation_->getNumberOfCells();
    int dimensionality = triangulation_->getDimensionality();
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(dimensionality + 1);

    for(SimplexId i = 0; i < cellNumber; i++) {
      for(int j = 0; j <= dimensionality; j++) {
        SimplexId vertexId = -1;
        triangulation_->getCellVertex(i, j, vertexId);
        idList->SetId(j, vertexId);
      }
      if(dimensionality == 2) {
        vtkUnstructuredGrid_->InsertNextCell(VTK_TRIANGLE, idList);
      } else if(dimensionality == 3) {
        vtkUnstructuredGrid_->InsertNextCell(VTK_TETRA, idList);
      }
    }

    {
      stringstream msg;
      msg << "[ttkTriangulation] Explicit triangulation computed in "
          << t.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  } else if((inputDataSet_->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)) {
    return vtkUnstructuredGrid::SafeDownCast(inputDataSet_);
  } else {
    stringstream msg;
    msg << "[ttkTriangulation] Unsupported input VTK class `"
        << inputDataSet_->GetClassName()
        << "' (ref=" << inputDataSet_->GetDataObjectType() << ")" << endl;
    msg << "[ttkTriangulation] Leaving an empty triangulation..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
  }

  return vtkUnstructuredGrid_.GetPointer();
}

bool ttkTriangulation::hasChangedConnectivity(Triangulation *triangulation,
                                              vtkDataSet *dataSet,
                                              vtkObject *callingObject) {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation)
    return false;
  if(!dataSet)
    return false;
  if(!callingObject)
    return false;
#endif

  if((dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)) {
    return (((vtkUnstructuredGrid *)dataSet)->GetCells()->GetMTime()
            > callingObject->GetMTime());
  } else if((dataSet->GetDataObjectType() == VTK_POLY_DATA)) {
    return (((vtkPolyData *)dataSet)->GetPolys()->GetMTime()
            > callingObject->GetMTime());
  } else if((dataSet->GetDataObjectType() == VTK_IMAGE_DATA)) {

    int vtkDimensions[3];
    vector<int> ttkDimensions;

    ((vtkImageData *)dataSet)->GetDimensions(vtkDimensions);
    triangulation->getGridDimensions(ttkDimensions);

    return !((vtkDimensions[0] == ttkDimensions[0])
             && (vtkDimensions[1] == ttkDimensions[1])
             && (vtkDimensions[2] == ttkDimensions[2]));
  } else {
    stringstream msg;
    Debug d;
    msg << "[ttkTriangulation] Unsupported input VTK class `"
        << dataSet->GetClassName() << "' (ref=" << dataSet->GetDataObjectType()
        << ")" << endl;
    d.dMsg(cerr, msg.str(), Debug::fatalMsg);
  }

  return true;
}

int ttkTriangulation::setInputData(vtkDataSet *dataSet) {

  int ret = 0;
  if(!triangulation_) {
    allocate();
  }

  if(dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID) {

    auto vtuDataSet = vtkUnstructuredGrid::SafeDownCast(dataSet);

    // Points
    if(vtuDataSet->GetPoints()) {
      switch(vtuDataSet->GetPoints()->GetDataType()) {
        case VTK_FLOAT:
          triangulation_->setInputPoints(
            vtuDataSet->GetNumberOfPoints(),
            ttkUtils::GetVoidPointer(vtuDataSet->GetPoints()->GetData()));
          break;
        case VTK_DOUBLE:
          triangulation_->setInputPoints(
            vtuDataSet->GetNumberOfPoints(),
            ttkUtils::GetVoidPointer(vtuDataSet->GetPoints()->GetData()), true);
          break;
        default:
          stringstream msg;
          msg << "[ttkTriangulation] Unsupported precision for input points!"
              << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          break;
      }
    }

    // Cells
    if(vtuDataSet->GetNumberOfCells()) {
      // As the documentation of this method says, we expect a triangulation
      // (simplicial complex in kD)
      std::string errMsg;

      // NOTE A mix of different simplices should theoretically be supported
      // throughout TTK but this appears to not hold up (yet) in practice.
      // Some more notes: https://github.com/topology-tool-kit/ttk/pull/323
      // maybe in the future instead of the current checks:
      //      ret = vtkSmartPointer<CellChecker>::New()->check3d(vtuDataSet);

      if(vtuDataSet->IsHomogeneous()) {
        const auto type = vtuDataSet->GetCellType(0);
        if(type != VTK_TETRA && type != VTK_TRIANGLE && type != VTK_LINE
           && type != VTK_VERTEX) {
          errMsg = "Input explicit mesh is not simplicial (at all). ";
          ret = -1;
        }
      } else {
        errMsg = "Input explicit mesh has different cell types. ";
        ret = -1;
      }

      if(ret != 0) {
        // NOTE In the long run it would be nice to not only detect this but
        // also fix it, for an improved user-experience. Maybe not at this place
        // though (or adapt doc) but earlier, like the macro
        // TTK_UNSTRUCTURED_GRID_NEW in ttkWrapper.h?

        // This class does not inherit from vtkObject, so the vtkWarningMacro
        // does not work. But no problem: `this` is only used to call
        // GetClassName which is not a virtual function so it would just return
        // "vtkObject" anyway.
        errMsg += "Please apply the Tetrahedralize filter before proceeding "
                  "with TTK. "
                  "The application behavior may be unspecified otherwise.";
        vtkWarningWithObjectMacro(nullptr, << errMsg);
      }

#ifdef TTK_CELL_ARRAY_NEW
      auto cellsConnectivity
        = static_cast<LongSimplexId *>(ttkUtils::GetVoidPointer(
          vtuDataSet->GetCells()->GetConnectivityArray(), 0));
      auto cellsOffsets = static_cast<LongSimplexId *>(
        ttkUtils::GetVoidPointer(vtuDataSet->GetCells()->GetOffsetsArray(), 0));
#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
      triangulation_->setInputCells(
        vtuDataSet->GetNumberOfCells(), cellsConnectivity, cellsOffsets);
#else
      auto extra_connectivity
        = reinterpret_cast<long long *>(cellsConnectivity);
      auto extra_offsets = reinterpret_cast<long long *>(cellsOffsets);
      triangulation_->setInputCells(
        dataSet->GetNumberOfCells(), extra_connectivity, extra_offsets);
#endif
#else
      auto cellsData = static_cast<LongSimplexId *>(
        ttkUtils::GetVoidPointer(vtuDataSet->GetCells()->GetData(), 0));
#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
      triangulation_->setInputCells(vtuDataSet->GetNumberOfCells(), cellsData);
#else
      auto extra_pt = reinterpret_cast<long long *>(cellsData);
      triangulation_->setInputCells(dataSet->GetNumberOfCells(), extra_pt);
#endif
#endif
    }
    inputDataSet_ = vtuDataSet;

  } else if(dataSet->GetDataObjectType() == VTK_POLY_DATA) {

    auto vtpDataSet = vtkPolyData::SafeDownCast(dataSet);

    // Points
    if(vtpDataSet->GetPoints()) {
      switch(vtpDataSet->GetPoints()->GetDataType()) {
        case VTK_FLOAT:
          triangulation_->setInputPoints(
            vtpDataSet->GetNumberOfPoints(),
            ttkUtils::GetVoidPointer(vtpDataSet->GetPoints()->GetData()));
          break;
        case VTK_DOUBLE:
          triangulation_->setInputPoints(
            vtpDataSet->GetNumberOfPoints(),
            ttkUtils::GetVoidPointer(vtpDataSet->GetPoints()->GetData()), true);
          break;
        default:
          stringstream msg;
          msg << "[ttkTriangulation] Unsupported precision for input points!"
              << endl;
          dMsg(cerr, msg.str(), Debug::fatalMsg);
          break;
      }
    }

    // Cells
    // NOTE Currently we do not handle vtkPolyData instances whose cells include
    // unconnected vertices or triangle strips
    //    if(polyData->GetNumberOfVerts())
    //      ; // add some 0D simplices to this ttkTriangulation
    //    if(polyData->GetNumberOfStrips())
    //      ; // add more 2D simplices to this ttkTriangulation (de-strip)
    if(vtpDataSet->GetNumberOfPolys() || vtpDataSet->GetNumberOfLines()) {
#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
      if(vtpDataSet->GetNumberOfPolys() > 0) { // have 2D

        ret = vtkSmartPointer<CellChecker>::New()->check2d(vtpDataSet);
        // NOTE If there are any 1D cells, we just ignore them
#ifdef TTK_CELL_ARRAY_NEW
        auto polyDataConnect
          = static_cast<LongSimplexId *>(ttkUtils::GetVoidPointer(
            vtpDataSet->GetPolys()->GetConnectivityArray(), 0));
        auto polyDataOffset
          = static_cast<LongSimplexId *>(ttkUtils::GetVoidPointer(
            vtpDataSet->GetPolys()->GetOffsetsArray(), 0));
        triangulation_->setInputCells(
          vtpDataSet->GetNumberOfCells(), polyDataConnect, polyDataOffset);
#else
        auto polysData = static_cast<LongSimplexId *>(
          ttkUtils::GetVoidPointer(vtpDataSet->GetPolys()->GetData(), 0));
        triangulation_->setInputCells(
          vtpDataSet->GetNumberOfCells(), polysData);
#endif

      } else if(vtpDataSet->GetNumberOfLines() > 0) { // have 1D

        ret = vtkSmartPointer<CellChecker>::New()->check1d(vtpDataSet);
#ifdef TTK_CELL_ARRAY_NEW
        auto lineDataConnect
          = static_cast<LongSimplexId *>(ttkUtils::GetVoidPointer(
            vtpDataSet->GetLines()->GetConnectivityArray(), 0));
        auto lineDataOffset
          = static_cast<LongSimplexId *>(ttkUtils::GetVoidPointer(
            vtpDataSet->GetLines()->GetOffsetsArray(), 0));
        triangulation_->setInputCells(
          vtpDataSet->GetNumberOfCells(), lineDataConnect, lineDataOffset);
#else
        auto linesData = static_cast<LongSimplexId *>(
          ttkUtils::GetVoidPointer(vtpDataSet->GetLines()->GetData(), 0));
        triangulation_->setInputCells(
          vtpDataSet->GetNumberOfCells(), linesData);
#endif
      } else {
        stringstream msg;
        msg << "[ttkTriangulation] empty polydata." << endl;
        dMsg(cerr, msg.str(), Debug::fatalMsg);
      }
#else
      int *pt = polyData->GetPolys()->GetPointer();
      if(!pt) {
        // 1D
        pt = polyData->GetLines()->GetPointer();
      }
      auto extra_pt
        = reinterpret_cast<long long *>(polysData ? polysData : linesData);
      triangulation_->setInputCells(dataSet->GetNumberOfCells(), extra_pt);
#endif
    }
    inputDataSet_ = vtpDataSet;

  } else if((dataSet->GetDataObjectType() == VTK_IMAGE_DATA)) {

    auto vtiData = vtkImageData::SafeDownCast(dataSet);

    int extents[6];
    vtiData->GetExtent(extents);

    double origin[3];
    vtiData->GetOrigin(origin);

    double spacing[3];
    vtiData->GetSpacing(spacing);

    int gridDimensions[3];
    vtiData->GetDimensions(gridDimensions);

    double firstPoint[3];
    firstPoint[0] = origin[0] + extents[0] * spacing[0];
    firstPoint[1] = origin[1] + extents[2] * spacing[1];
    firstPoint[2] = origin[2] + extents[4] * spacing[2];

    triangulation_->setInputGrid(
      firstPoint[0], firstPoint[1], firstPoint[2], spacing[0], spacing[1],
      spacing[2], gridDimensions[0], gridDimensions[1], gridDimensions[2]);

    inputDataSet_ = dataSet;
  } else {
    stringstream msg;
    msg << "[ttkTriangulation] Unsupported input VTK class `"
        << dataSet->GetClassName() << "' (ref=" << dataSet->GetDataObjectType()
        << ")" << endl;
    msg << "[ttkTriangulation] Leaving an empty triangulation..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
  }

  return ret;
}

int ttkTriangulation::shallowCopy(vtkDataObject *other) {

  if((triangulation_) && (hasAllocated_)) {
    delete triangulation_;
  }
  triangulation_ = NULL;

  if((other) && (((vtkDataSet *)other)->GetNumberOfPoints())) {

    ttkTriangulation *otherTriangulation
      = dynamic_cast<ttkTriangulation *>(other);

    if(otherTriangulation) {
      triangulation_ = otherTriangulation->triangulation_;
      hasAllocated_ = false;
    } else {
      // let's create the object
      allocate();
    }

    // populate the data-structure
    if((triangulation_) && (triangulation_->isEmpty())) {
      setInputData((vtkDataSet *)other);
    }
  }

  return 0;
}

vtkStandardNewMacro(ttkUnstructuredGrid)

  ttkUnstructuredGrid::ttkUnstructuredGrid() {
}

ttkUnstructuredGrid::~ttkUnstructuredGrid() {
}

void ttkUnstructuredGrid::CopyStructure(vtkDataSet *other) {

  vtkUnstructuredGrid::CopyStructure(other);
  ttkTriangulation::shallowCopy(other);
}

void ttkUnstructuredGrid::DeepCopy(vtkDataObject *other) {

  vtkUnstructuredGrid::DeepCopy(other);
  ttkTriangulation::deepCopy(other);
}

void ttkUnstructuredGrid::ShallowCopy(vtkDataObject *other) {

  vtkUnstructuredGrid::ShallowCopy(other);
  ttkTriangulation::shallowCopy(other);
}

vtkStandardNewMacro(ttkImageData)

  ttkImageData::ttkImageData() {
}

ttkImageData::~ttkImageData() {
}

void ttkImageData::CopyStructure(vtkDataSet *other) {

  vtkImageData::CopyStructure(other);
  ttkTriangulation::shallowCopy(other);
}

void ttkImageData::DeepCopy(vtkDataObject *other) {

  vtkImageData::DeepCopy(other);
  ttkTriangulation::deepCopy(other);
}

void ttkImageData::ShallowCopy(vtkDataObject *other) {

  vtkImageData::ShallowCopy(other);
  ttkTriangulation::shallowCopy(other);
}

vtkStandardNewMacro(ttkPolyData)

  ttkPolyData::ttkPolyData() {
}

ttkPolyData::~ttkPolyData() {
}

void ttkPolyData::CopyStructure(vtkDataSet *other) {

  vtkPolyData::CopyStructure(other);
  ttkTriangulation::shallowCopy(other);
}

void ttkPolyData::DeepCopy(vtkDataObject *other) {

  vtkPolyData::DeepCopy(other);
  ttkTriangulation::deepCopy(other);
}

void ttkPolyData::ShallowCopy(vtkDataObject *other) {

  vtkPolyData::ShallowCopy(other);
  ttkTriangulation::shallowCopy(other);
}
