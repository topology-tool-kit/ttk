#include <ttkTriangulationFactory.h>

#include <Triangulation.h>
#include <ttkUtils.h>
#include <vtkCellTypes.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

vtkCellArray *GetCells(vtkDataSet *dataSet) {
  switch(dataSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID: {
      auto dataSetAsUG = static_cast<vtkUnstructuredGrid *>(dataSet);
      return dataSetAsUG->GetCells();
    }
    case VTK_POLY_DATA: {
      auto dataSetAsPD = static_cast<vtkPolyData *>(dataSet);
      return dataSetAsPD->GetNumberOfPolys() > 0   ? dataSetAsPD->GetPolys()
             : dataSetAsPD->GetNumberOfLines() > 0 ? dataSetAsPD->GetLines()
                                                   : dataSetAsPD->GetVerts();
    }
  }
  return nullptr;
}

int checkCellTypes(vtkPointSet *object) {
  auto cellTypes = vtkSmartPointer<vtkCellTypes>::New();
  object->GetCellTypes(cellTypes);

  size_t nTypes = cellTypes->GetNumberOfTypes();

  // if cells are empty
  if(nTypes == 0)
    return 1; // no error

  // if cells are not homogeneous
  if(nTypes > 1)
    return -1;

  // if cells are not simplices
  if(nTypes == 1) {
    const auto &cellType = cellTypes->GetCellType(0);
    if(cellType != VTK_VERTEX && cellType != VTK_LINE
       && cellType != VTK_TRIANGLE && cellType != VTK_TETRA)
      return -2;
  }

  return 1;
}

struct ttkOnDeleteCommand : public vtkCommand {
  RegistryKey key;
  vtkObject *observee;

  static ttkOnDeleteCommand *New() {
    return new ttkOnDeleteCommand;
  }
  vtkTypeMacro(ttkOnDeleteCommand, vtkCommand);

  void Init(vtkDataSet *dataSet) {
    this->key = ttkTriangulationFactory::GetKey(dataSet);

    if(dataSet->IsA("vtkPointSet"))
      this->observee = static_cast<vtkObject *>(GetCells(dataSet));
    else
      this->observee = static_cast<vtkObject *>(dataSet);

    this->observee->AddObserver(vtkCommand::DeleteEvent, this, 1);
  }

  void Execute(vtkObject *,
               unsigned long ttkNotUsed(eventId),
               void *ttkNotUsed(callData)) override {
    if(this->observee)
      this->observee->RemoveObserver(this);

    auto instance = &ttkTriangulationFactory::Instance;

    if(instance->registry.empty()) {
      return;
    }

    auto it = instance->registry.find(this->key);
    if(it != instance->registry.end()) {
      instance->registry.erase(it);
      instance->printMsg("Triangulation Deleted", ttk::debug::Priority::DETAIL);
      instance->printMsg("# Registered Triangulations: "
                           + std::to_string(instance->registry.size()),
                         ttk::debug::Priority::VERBOSE);
    }
  }
};

RegistryValue::RegistryValue(vtkDataSet *dataSet,
                             ttk::Triangulation *triangulation_)
  : triangulation(triangulation_), owner(dataSet) {
  auto cells = GetCells(dataSet);
  if(cells)
    this->cellModTime = cells->GetMTime();

  if(dataSet->IsA("vtkImageData")) {
    auto image = static_cast<vtkImageData *>(dataSet);
    image->GetExtent(this->extent);
    image->GetOrigin(this->origin);
    image->GetSpacing(this->spacing);
    image->GetDimensions(this->dimensions);
  }

  auto onDelete = vtkSmartPointer<ttkOnDeleteCommand>::New();
  onDelete->Init(dataSet);
}

bool RegistryValue::isValid(vtkDataSet *dataSet) const {
  auto cells = GetCells(dataSet);
  if(cells)
    return this->cellModTime == cells->GetMTime();

  if(dataSet->IsA("vtkImageData")) {
    auto image = static_cast<vtkImageData *>(dataSet);

    int extent_[6];
    double origin_[3];
    double spacing_[3];
    int dimensions_[3];

    image->GetExtent(extent_);
    image->GetOrigin(origin_);
    image->GetSpacing(spacing_);
    image->GetDimensions(dimensions_);

    bool isValid = true;
    for(int i = 0; i < 6; i++)
      if(this->extent[i] != extent_[i])
        isValid = false;
    for(int i = 0; i < 3; i++)
      if(this->origin[i] != origin_[i] || this->spacing[i] != spacing_[i]
         || this->dimensions[i] != dimensions_[i])
        isValid = false;

    return isValid;
  }

  return false;
}

ttkTriangulationFactory::ttkTriangulationFactory() {
  this->setDebugMsgPrefix("TriangulationFactory");
}

RegistryTriangulation
  ttkTriangulationFactory::CreateImplicitTriangulation(vtkImageData *image) {
  ttk::Timer timer;
  this->printMsg("Initializing Implicit Triangulation", 0, 0,
                 ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

  auto triangulation = RegistryTriangulation(new ttk::Triangulation());

  int extent[6];
  image->GetExtent(extent);

  double origin[3];
  image->GetOrigin(origin);

  double spacing[3];
  image->GetSpacing(spacing);

  int dimensions[3];
  image->GetDimensions(dimensions);

  double firstPoint[3];
  firstPoint[0] = origin[0] + extent[0] * spacing[0];
  firstPoint[1] = origin[1] + extent[2] * spacing[1];
  firstPoint[2] = origin[2] + extent[4] * spacing[2];

  triangulation->setInputGrid(firstPoint[0], firstPoint[1], firstPoint[2],
                              spacing[0], spacing[1], spacing[2], dimensions[0],
                              dimensions[1], dimensions[2]);

  this->printMsg("Initializing Implicit Triangulation", 1,
                 timer.getElapsedTime(), ttk::debug::LineMode::NEW,
                 ttk::debug::Priority::DETAIL);

  return triangulation;
}

RegistryTriangulation
  ttkTriangulationFactory::CreateExplicitTriangulation(vtkPointSet *pointSet) {
  ttk::Timer timer;

  auto points = pointSet->GetPoints();
  if(!points) {
    this->printErr("DataSet has uninitialized `vtkPoints`.");
    return nullptr;
  }

  auto cells = GetCells(pointSet);
  if(!cells) {
    this->printErr("DataSet has uninitialized `vtkCellArray`.");
    return nullptr;
  }

  auto triangulation = RegistryTriangulation(new ttk::Triangulation());
  int hasIndexArray
    = pointSet->GetPointData()->HasArray(ttk::compactTriangulationIndex);

  if(hasIndexArray) {
    this->printMsg("Initializing Compact Triangulation", 0, 0,
                   ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);
  } else {
    this->printMsg("Initializing Explicit Triangulation", 0, 0,
                   ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);
  }

  // Points
  {
    auto pointDataType = points->GetDataType();
    if(pointDataType != VTK_FLOAT && pointDataType != VTK_DOUBLE) {
      this->printErr("Unable to initialize 'ttk::Triangulation' for point "
                     "precision other than 'float' or 'double'.");
      return {};
    }

    void *pointDataArray = ttkUtils::GetVoidPointer(points);
    if(hasIndexArray) {
      vtkAbstractArray *indexArray = pointSet->GetPointData()->GetAbstractArray(
        ttk::compactTriangulationIndex);
      triangulation->setStellarInputPoints(
        points->GetNumberOfPoints(), pointDataArray,
        (int *)indexArray->GetVoidPointer(0), pointDataType == VTK_DOUBLE);
    } else {
      triangulation->setInputPoints(points->GetNumberOfPoints(), pointDataArray,
                                    pointDataType == VTK_DOUBLE);
    }
  }

  // check if cell types are simplices
  int cellTypeStatus = checkCellTypes(pointSet);
  if(cellTypeStatus == -1) {
    this->printWrn("Inhomogeneous cell dimensions detected.");
    this->printWrn(
      "Consider using `ttkExtract` to extract cells of a given dimension.");
    return {};
  } else if(cellTypeStatus == -2) {
    this->printWrn("Cells are not simplices.");
    this->printWrn("Consider using `vtkTetrahedralize` in pre-processing.");
    return {};
  }

  // Cells
  int nCells = cells->GetNumberOfCells();
  if(nCells > 0) {
    if(!cells->IsStorage64Bit()) {
      if(cells->CanConvertTo64BitStorage()) {
        this->printWrn("Converting the cell array to 64-bit storage");
        bool success = cells->ConvertTo64BitStorage();
        if(!success) {
          this->printErr(
            "Error converting the provided cell array to 64-bit storage");
          return {};
        }
      } else {
        this->printErr(
          "Cannot convert the provided cell array to 64-bit storage");
        return {};
      }
    }
    auto connectivity = static_cast<vtkIdType *>(
      ttkUtils::GetVoidPointer(cells->GetConnectivityArray()));
    auto offsets = static_cast<vtkIdType *>(
      ttkUtils::GetVoidPointer(cells->GetOffsetsArray()));

    int status;
    if(hasIndexArray) {
      status
        = triangulation->setStellarInputCells(nCells, connectivity, offsets);
    } else {
      status = triangulation->setInputCells(nCells, connectivity, offsets);
    }

    if(status != 0) {
      this->printErr(
        "Run the `vtkTetrahedralize` filter to resolve the issue.");
      return {};
    }
  }

  if(hasIndexArray) {
    this->printMsg("Initializing Compact Triangulation", 1,
                   timer.getElapsedTime(), ttk::debug::LineMode::NEW,
                   ttk::debug::Priority::DETAIL);
  } else {
    this->printMsg("Initializing Explicit Triangulation", 1,
                   timer.getElapsedTime(), ttk::debug::LineMode::NEW,
                   ttk::debug::Priority::DETAIL);
  }

  return triangulation;
}

RegistryTriangulation
  ttkTriangulationFactory::CreateTriangulation(vtkDataSet *dataSet) {
  switch(dataSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
    case VTK_POLY_DATA: {
      return this->CreateExplicitTriangulation(
        static_cast<vtkPointSet *>(dataSet));
    }
    case VTK_IMAGE_DATA: {
      return this->CreateImplicitTriangulation((vtkImageData *)dataSet);
    }
    default: {
      this->printErr("Unable to triangulate `"
                     + std::string(dataSet->GetClassName()) + "`");
    }
  }

  return nullptr;
}

ttk::Triangulation *ttkTriangulationFactory::GetTriangulation(
  int debugLevel, float cacheRatio, vtkDataSet *object) {
  auto instance = &ttkTriangulationFactory::Instance;
  instance->setDebugLevel(debugLevel);

  auto key = ttkTriangulationFactory::GetKey(object);

  ttk::Triangulation *triangulation{nullptr};
  auto it = instance->registry.find(key);
  if(it != instance->registry.end()) {
    // object is the owner of the explicit or implicit triangulation
    if(it->second.isValid(object)) {
      instance->printMsg(
        "Retrieving Existing Triangulation", ttk::debug::Priority::DETAIL);
      triangulation = it->second.triangulation.get();
    } else {
      instance->printMsg(
        "Existing Triangulation No Longer Valid", ttk::debug::Priority::DETAIL);
      instance->registry.erase(key);
    }
  }

  if(!triangulation && object->IsA("vtkImageData")) {
    instance->FindImplicitTriangulation(
      triangulation, static_cast<vtkImageData *>(object));
    if(triangulation)
      instance->printMsg("Retrieving Equivalent Implicit-Triangulation",
                         ttk::debug::Priority::DETAIL);
  }

  if(!triangulation) {
    triangulation = instance->CreateTriangulation(object).release();
    if(triangulation) {
      instance->registry.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(key),
                                 std::forward_as_tuple(object, triangulation));
    }
  }

  instance->printMsg(
    "# Registered Triangulations: " + std::to_string(instance->registry.size()),
    ttk::debug::Priority::VERBOSE);

  if(triangulation) {
    triangulation->setDebugLevel(debugLevel);
    triangulation->setCacheSize(cacheRatio);
  }

  return triangulation;
}

int ttkTriangulationFactory::FindImplicitTriangulation(
  ttk::Triangulation *&triangulation, vtkImageData *image) {

  for(const auto &it : this->registry) {
    if(it.second.owner->IsA("vtkImageData")) {
      if(it.second.isValid(image)) {
        triangulation = it.second.triangulation.get();
        return 1;
      }
    }
  }

  return 0;
}

RegistryKey ttkTriangulationFactory::GetKey(vtkDataSet *dataSet) {
  switch(dataSet->GetDataObjectType()) {
    case VTK_IMAGE_DATA: {
      return (RegistryKey)dataSet;
    }
    default: {
      auto cells = GetCells(dataSet);
      if(cells)
        return (RegistryKey)cells;
    }
  }

  return 0;
}

ttkTriangulationFactory ttkTriangulationFactory::Instance{};
