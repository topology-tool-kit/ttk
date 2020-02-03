#include <ttkTriangulation.h>
#include <ttkTriangulationAlgorithm.h>

using namespace std;
using namespace ttk;

/**
 * @brief Auxiliary struct to check if VTK cells are simplices.
 */
struct CellChecker : public vtkObject // inherit to be able to use vtkErrorMacro
{
  vtkTypeMacro(CellChecker, vtkObject)

  CellChecker() = default;
  static CellChecker* New() { return new CellChecker; }

  static std::string getErrMsg() // in C++ > 11: static constexpr errMsg = "foo";
  {
    return "Input explicit mesh is not (completely) simplicial. "
           "Please apply the Tetrahedralize filter before any TTK filter. "
           "The application may now crash anytime or produce false results.";
  }

  int check1d(vtkDataSet* ds)
  {
    // e.g. VTK_POLY_LINE is not supported (as from vtkPolyData::GetLines)
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c)
    {
      if(ds->GetCellType(c) != VTK_LINE)
      {
        vtkErrorMacro(<< getErrMsg())
        return -1;
      }
    }
    return 0;
  }

  int check2d(vtkDataSet* ds)
  {
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c)
    {
      const auto type = ds->GetCellType(c);
      if(type != VTK_TRIANGLE && type != VTK_LINE)
      {
        vtkErrorMacro(<< getErrMsg())
        return -1;
      }
    }
    return 0;
  }

  int check3d(vtkDataSet* ds)
  {
    const vtkIdType nc = ds->GetNumberOfCells();
    for(vtkIdType c = 0; c < nc; ++c)
    {
      const auto type = ds->GetCellType(c);
      if(type != VTK_TETRA && type != VTK_TRIANGLE && type != VTK_LINE)
      {
        vtkErrorMacro(<< getErrMsg())
        return -1;
      }
    }
    return 0;
  }
};

ttkTriangulation::ttkTriangulation() {

  inputDataSet_ = NULL;

  hasAllocated_ = false;
  triangulation_ = NULL;
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

  if((dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)) {

    const auto grid = static_cast<vtkUnstructuredGrid *>(dataSet);
    if(grid->GetPoints()) {
      if(grid->GetPoints()->GetDataType() == VTK_FLOAT) {
        triangulation_->setInputPoints(
          dataSet->GetNumberOfPoints(), grid->GetPoints()->GetVoidPointer(0));
      } else if(grid->GetPoints()->GetDataType() == VTK_DOUBLE) {
        triangulation_->setInputPoints(dataSet->GetNumberOfPoints(),
                                       grid->GetPoints()->GetVoidPointer(0),
                                       true);
      } else {
        stringstream msg;
        msg << "[ttkTriangulation] Unsupported precision for input points!"
            << endl;
        dMsg(cerr, msg.str(), Debug::fatalMsg);
      }
    }
    if(grid->GetCells()) {
      // As the documentation of this method says, we expect a triangulation
      // (simplicial complex in kD)
      auto checkSimplicial = [](int type){
        return type == VTK_LINE || type == VTK_TRIANGLE || type == VTK_TETRA;
      };
      std::string errMsg;

      // NOTE A mix of different simplices should theoretically be supported
      // but this appears to not hold up (yet) in practice. Current discussion:
      // https://github.com/topology-tool-kit/ttk/pull/323
//      ret = vtkSmartPointer<CellChecker>::New()->check3d(grid); // maybe in the future

      if(grid->IsHomogeneous()) {
        const auto type = grid->GetCellType(0);
        if(type != VTK_TETRA && type != VTK_TRIANGLE && type != VTK_LINE) {
          errMsg = "Input explicit mesh is not simplicial (at all). ";
          ret = -1;
        }
      } else
      {
        errMsg = "Input explicit mesh has different cell types. ";
        ret = -1;
      }

      if(ret != 0) {
        // NOTE In the long run it would be nice to not only detect this but also
        // fix it, for an improved user-experience.
        // Maybe not at this place though (or adapt doc) but earlier like the macro
        // TTK_UNSTRUCTURED_GRID_NEW in ttkWrapper.h?

        // This class does not inherit from vtkObject, so the vtkErrorMacro does not work.
        // But no problem: `this` is only used to call GetClassName which is not a virtual
        // function so it would just return "vtkObject" anyway.
        errMsg += "Please apply the Tetrahedralize filter before any TTK filter. "
                  "The application may now crash anytime or produce false results.";
        vtkErrorWithObjectMacro(nullptr, << errMsg)
      }

#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
      triangulation_->setInputCells(
        dataSet->GetNumberOfCells(), grid->GetCells()->GetPointer());
#else
      int *pt = grid->GetCells()->GetPointer();
      long long extra_pt = *pt;
      triangulation_->setInputCells(dataSet->GetNumberOfCells(), &extra_pt);
#endif
    }
    inputDataSet_ = dataSet;
  } else if((dataSet->GetDataObjectType() == VTK_POLY_DATA)) {

    const auto polyData = static_cast<vtkPolyData*>(dataSet);
    if(polyData->GetPoints()) {
      if(polyData->GetPoints()->GetDataType() == VTK_FLOAT) {
        triangulation_->setInputPoints(
          dataSet->GetNumberOfPoints(), polyData->GetPoints()->GetVoidPointer(0));
      } else if(polyData->GetPoints()->GetDataType() == VTK_DOUBLE) {
        triangulation_->setInputPoints(
          dataSet->GetNumberOfPoints(), polyData->GetPoints()->GetVoidPointer(0), true);
      } else {
        stringstream msg;
        msg << "[ttkTriangulation] Unsupported precision for input points!"
            << endl;
        dMsg(cerr, msg.str(), Debug::fatalMsg);
      }
      triangulation_->setInputPoints(
        dataSet->GetNumberOfPoints(), (float *)polyData->GetPoints()->GetVoidPointer(0));
    }

    // NOTE Currently we do not handle vtkPolyData instances whose cells include
    // unconnected vertices or triangle strips
//    if(polyData->GetNumberOfVerts())
//      ;
//    if(polyData->GetNumberOfStrips())
//      ;

    if(polyData->GetPolys()) {
#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
      // does never return a null pointer, only an empty array
//      if(polyData->GetPolys()->GetPointer()) {
      if(polyData->GetNumberOfPolys()) { // have 2D
        ret = vtkSmartPointer<CellChecker>::New()->check2d(polyData);
        // NOTE If there are any 1D cells, we just ignore them
        triangulation_->setInputCells(
          dataSet->GetNumberOfCells(), polyData->GetPolys()->GetPointer());
      } else if(polyData->GetNumberOfLines()) { // have 1D
        ret = vtkSmartPointer<CellChecker>::New()->check1d(polyData);
        triangulation_->setInputCells(
          dataSet->GetNumberOfCells(), polyData->GetLines()->GetPointer());
      }
#else
      int *pt = polyData->GetPolys()->GetPointer();
      if(!pt) {
        // 1D
        pt = polyData->GetLines()->GetPointer();
      }
      long long extra_pt = *pt;
      triangulation_->setInputCells(dataSet->GetNumberOfCells(), &extra_pt);
#endif
    }
    inputDataSet_ = dataSet;
  } else if((dataSet->GetDataObjectType() == VTK_IMAGE_DATA)) {
    vtkImageData *imageData = (vtkImageData *)dataSet;

    int extents[6];
    imageData->GetExtent(extents);

    double origin[3];
    imageData->GetOrigin(origin);

    double spacing[3];
    imageData->GetSpacing(spacing);

    int gridDimensions[3];
    imageData->GetDimensions(gridDimensions);

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
