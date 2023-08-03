#include <ttkMacros.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h> // for VTK_VERSION_CHECK via ParaView 5.8.1

int VTUToDiagram(ttk::DiagramType &diagram,
                 vtkUnstructuredGrid *vtu,
                 const ttk::Debug &dbg) {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();
  const auto points = vtu->GetPoints();

  if(pd == nullptr) {
    dbg.printErr("VTU diagram with NULL Point Data");
    return -1;
  }
  if(cd == nullptr) {
    dbg.printErr("VTU diagram with NULL Cell Data");
    return -2;
  }
  if(points == nullptr) {
    dbg.printErr("VTU with no points");
    return -3;
  }

  // cell data
  const auto pairId = vtkIntArray::SafeDownCast(
    cd->GetArray(ttk::PersistencePairIdentifierName));
  const auto pairType
    = vtkIntArray::SafeDownCast(cd->GetArray(ttk::PersistencePairTypeName));
  const auto pairPers = cd->GetArray(ttk::PersistenceName);
  const auto birthScalars = cd->GetArray(ttk::PersistenceBirthName);
  const auto isFinite = cd->GetArray(ttk::PersistenceIsFinite);

  // point data
  const auto vertexId
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto critType
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::PersistenceCriticalTypeName));
  const auto coords = vtkFloatArray::SafeDownCast(
    pd->GetArray(ttk::PersistenceCoordinatesName));

  const bool embed = coords == nullptr;

  if(pairId == nullptr) {
    dbg.printErr("Missing PairIdentifier cell data array");
    return -5;
  }
  if(pairType == nullptr) {
    dbg.printErr("Missing PairType cell data array");
    return -6;
  }
  if(pairPers == nullptr) {
    dbg.printErr("Missing Persistence cell data array");
    return -7;
  }
  if(vertexId == nullptr) {
    dbg.printErr("Missing ttkVertexScalarField point data array");
    return -8;
  }
  if(critType == nullptr) {
    dbg.printErr("Missing CriticalType point data array");
    return -9;
  }
  if(birthScalars == nullptr) {
    dbg.printErr("Missing Birth cell data array");
    return -10;
  }
  if(isFinite == nullptr) {
    dbg.printErr("Missing IsFinite cell data array");
    return -12;
  }

  int nPairs = pairId->GetNumberOfTuples();

  // compact pairIds in [0, nPairs - 1] (diagonal excepted)
  for(int i = 0; i < nPairs; i++) {
    if(pairId->GetTuple1(i) != -1) {
      pairId->SetTuple1(i, i);
    } else {
      // detect diagram diagonal
      nPairs -= 1;
    }
  }

  if(nPairs < 1) {
    dbg.printErr("Diagram has no pairs");
    return -4;
  }

  diagram.resize(nPairs);

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < nPairs; ++i) {

    const auto pId = pairId->GetValue(i);
    // skip diagram diagonal if present
    if(pId == -1) {
      continue;
    }

    const auto i0 = 2 * i + 0;
    const auto i1 = 2 * i + 1;

    const auto v0 = vertexId->GetValue(i0);
    const auto v1 = vertexId->GetValue(i1);
    const auto ct0 = static_cast<ttk::CriticalType>(critType->GetValue(i0));
    const auto ct1 = static_cast<ttk::CriticalType>(critType->GetValue(i1));

    const auto pType = pairType->GetValue(i);
    const auto pers = pairPers->GetTuple1(i);
    const auto birth = birthScalars->GetTuple1(i);
    const auto isFin = static_cast<bool>(isFinite->GetTuple1(i));

    std::array<float, 3> coordsBirth{}, coordsDeath{};
    // no vtkPoints::GetPoint() taking a float array, have to do the
    // conversion by hand...
    std::array<double, 3> tmp{};

    if(embed) {
      points->GetPoint(i0, tmp.data());
      coordsBirth = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
      points->GetPoint(i1, tmp.data());
      coordsDeath = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
    } else {
      coords->GetTuple(i0, tmp.data());
      coordsBirth = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
      coords->GetTuple(i1, tmp.data());
      coordsDeath = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
    }

    // put pairs in diagram
    diagram[i] = ttk::PersistencePair{
      ttk::CriticalVertex{v0, ct0, birth, coordsBirth},
      ttk::CriticalVertex{v1, ct1, birth + pers, coordsDeath}, pType, isFin};
  }

  return 0;
}

int DiagramToVTU(vtkUnstructuredGrid *vtu,
                 const ttk::DiagramType &diagram,
                 vtkDataArray *const inputScalars,
                 const ttk::Debug &dbg,
                 const int dim,
                 const bool embedInDomain) {

  if(diagram.empty()) {
    dbg.printErr("Empty diagram");
    return -1;
  }

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  if(pd == nullptr || cd == nullptr) {
    dbg.printErr("Grid has no point data or no cell data");
    return -2;
  }

  // point data arrays

  vtkNew<ttkSimplexIdTypeArray> vertsId{};
  vertsId->SetName(ttk::VertexScalarFieldName);
  vertsId->SetNumberOfTuples(2 * diagram.size());
  pd->AddArray(vertsId);

  vtkNew<vtkIntArray> critType{};
  critType->SetName(ttk::PersistenceCriticalTypeName);
  critType->SetNumberOfTuples(2 * diagram.size());
  pd->AddArray(critType);

  vtkNew<vtkFloatArray> coordsScalars{};

  if(!embedInDomain) {
    coordsScalars->SetNumberOfComponents(3);
    coordsScalars->SetName(ttk::PersistenceCoordinatesName);
    coordsScalars->SetNumberOfTuples(2 * diagram.size());
    pd->AddArray(coordsScalars);
  }

  // cell data arrays

  vtkNew<ttkSimplexIdTypeArray> pairsId{};
  pairsId->SetName(ttk::PersistencePairIdentifierName);
  pairsId->SetNumberOfTuples(diagram.size());
  cd->AddArray(pairsId);

  vtkNew<vtkIntArray> pairsDim{};
  pairsDim->SetName(ttk::PersistencePairTypeName);
  pairsDim->SetNumberOfTuples(diagram.size());
  cd->AddArray(pairsDim);

  vtkSmartPointer<vtkDataArray> const persistence{inputScalars->NewInstance()};
  persistence->SetName(ttk::PersistenceName);
  persistence->SetNumberOfTuples(diagram.size());
  cd->AddArray(persistence);

  vtkSmartPointer<vtkDataArray> const birthScalars{inputScalars->NewInstance()};
  birthScalars->SetName(ttk::PersistenceBirthName);
  birthScalars->SetNumberOfTuples(diagram.size());
  cd->AddArray(birthScalars);

  vtkNew<vtkUnsignedCharArray> isFinite{};
  isFinite->SetName(ttk::PersistenceIsFinite);
  isFinite->SetNumberOfTuples(diagram.size());
  cd->AddArray(isFinite);

  // grid

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * diagram.size());
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(diagram.size() + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * diagram.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(dbg.getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &pair{diagram[i]};
    const auto i0{2 * i + 0}, i1{2 * i + 1};
    if(embedInDomain) {
      points->SetPoint(
        i0, pair.birth.coords[0], pair.birth.coords[1], pair.birth.coords[2]);
      points->SetPoint(
        i1, pair.death.coords[0], pair.death.coords[1], pair.death.coords[2]);
    } else {
      points->SetPoint(i0, pair.birth.sfValue, pair.birth.sfValue, 0);
      points->SetPoint(i1, pair.birth.sfValue, pair.death.sfValue, 0);
    }

    connectivity->SetTuple1(i0, i0);
    connectivity->SetTuple1(i1, i1);
    offsets->SetTuple1(i, 2 * i);

    // point data
    vertsId->SetTuple1(i0, pair.birth.id);
    vertsId->SetTuple1(i1, pair.death.id);
    critType->SetTuple1(i0, static_cast<ttk::SimplexId>(pair.birth.type));
    critType->SetTuple1(i1, static_cast<ttk::SimplexId>(pair.death.type));

    if(!embedInDomain) {
      coordsScalars->SetTuple3(
        i0, pair.birth.coords[0], pair.birth.coords[1], pair.birth.coords[2]);
      coordsScalars->SetTuple3(
        i1, pair.death.coords[0], pair.death.coords[1], pair.death.coords[2]);
    }

    // cell data
    pairsId->SetTuple1(i, i);
    persistence->SetTuple1(i, pair.persistence());
    birthScalars->SetTuple1(i, pair.birth.sfValue);
    isFinite->SetTuple1(i, pair.isFinite);
    pairsDim->SetTuple1(
      i, (pair.dim == 2 && pair.isFinite) ? dim - 1 : pair.dim);
  }
  offsets->SetTuple1(diagram.size(), connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  vtu->SetPoints(points);
  vtu->SetCells(VTK_LINE, cells);

  if(!embedInDomain) {
    // highest birth (last pair)
    const auto lastPair = std::max_element(diagram.begin(), diagram.end());
    // add diagonal (first point -> last birth/penultimate point)
    std::array<vtkIdType, 2> diag{
      0, 2 * std::distance(diagram.begin(), lastPair)};
    vtu->InsertNextCell(VTK_LINE, 2, diag.data());
    pairsId->InsertTuple1(diagram.size(), -1);
    pairsDim->InsertTuple1(diagram.size(), -1);
    isFinite->InsertTuple1(diagram.size(), false);
    // persistence of global min-max pair
    const auto maxPersistence = diagram[0].persistence();
    persistence->InsertTuple1(diagram.size(), 2 * maxPersistence);
    // birth == death == 0
    birthScalars->InsertTuple1(diagram.size(), 0);
  }

  return 0;
}

int ProjectDiagramInsideDomain(vtkUnstructuredGrid *const inputDiagram,
                               vtkUnstructuredGrid *const outputDiagram,
                               const ttk::Debug &dbg) {

  ttk::Timer tm{};

  // use vtkThreshold to remove diagonal (PairIdentifier == -1)
  vtkNew<vtkThreshold> threshold{};
  threshold->SetInputDataObject(0, inputDiagram);
  threshold->SetInputArrayToProcess(0, 0, 0,
                                    vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                    ttk::PersistencePairIdentifierName);
#if VTK_VERSION_NUMBER < VTK_VERSION_CHECK(9, 2, 0)
  threshold->ThresholdByUpper(0);
#else
  threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
  threshold->SetUpperThreshold(0);
#endif

  threshold->Update();

  auto diagonalLess = threshold->GetOutput();
  auto diagonalLessData = diagonalLess->GetPointData();

  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    diagonalLessData->GetAbstractArray(ttk::PersistenceCoordinatesName));

  // set new points from Coordinates array
  vtkNew<vtkFloatArray> coords{};
  coords->DeepCopy(critCoordinates);
  coords->SetName("Points");
  diagonalLess->GetPoints()->SetData(coords);
  diagonalLessData->RemoveArray(ttk::PersistenceCoordinatesName);

  outputDiagram->ShallowCopy(diagonalLess);

  // don't forget to forward the Field Data
  outputDiagram->GetFieldData()->ShallowCopy(inputDiagram->GetFieldData());

  dbg.printMsg("Projected Persistence Diagram inside domain", 1.0,
               tm.getElapsedTime(), dbg.getThreadNumber());

  return 0;
}

template <typename dataType>
void getCoords(vtkPoints *points,
               const dataType *const births,
               const dataType *const perss,
               vtkIdType nPoints,
               const int nThreads) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < nPoints / 2; ++i) {
    std::array<float, 3> pt0{
      static_cast<float>(births[i]),
      static_cast<float>(births[i]),
      0,
    };
    std::array<float, 3> pt1{
      static_cast<float>(births[i]),
      static_cast<float>(births[i] + perss[i]),
      0,
    };
    points->SetPoint(2 * i + 0, pt0.data());
    points->SetPoint(2 * i + 1, pt1.data());
  }

  TTK_FORCE_USE(nThreads);
}

int ProjectDiagramIn2D(vtkUnstructuredGrid *const inputDiagram,
                       vtkUnstructuredGrid *const outputDiagram,
                       const ttk::Debug &dbg) {

  ttk::Timer tm{};

  outputDiagram->ShallowCopy(inputDiagram);

  auto pointData = outputDiagram->GetPointData();

  auto birth = inputDiagram->GetCellData()->GetArray(ttk::PersistenceBirthName);
  auto pers = inputDiagram->GetCellData()->GetArray(ttk::PersistenceName);

  if(birth == nullptr || pers == nullptr) {
    dbg.printErr("Missing Birth or Persistence arrays");
    return 1;
  }

  // generate a new `Coordinates` pointData array
  vtkNew<vtkFloatArray> coords{};
  coords->DeepCopy(inputDiagram->GetPoints()->GetData());
  coords->SetName(ttk::PersistenceCoordinatesName);
  pointData->AddArray(coords);

  vtkNew<vtkPoints> points{};
  const auto nPoints = inputDiagram->GetNumberOfPoints();
  points->SetNumberOfPoints(nPoints);

  if(birth->GetNumberOfTuples() != nPoints / 2
     || pers->GetNumberOfTuples() != nPoints / 2) {
    dbg.printErr("Wrong number of tuples for Birth or Persistence arrays");
    return 2;
  }

  switch(birth->GetDataType()) {
    vtkTemplateMacro(getCoords(points, ttkUtils::GetPointer<VTK_TT>(birth),
                               ttkUtils::GetPointer<VTK_TT>(pers), nPoints,
                               dbg.getThreadNumber()));
  }

  outputDiagram->SetPoints(points);

  // add diagonal(first point -> last birth/penultimate point)
  std::array<vtkIdType, 2> diag{0, 2 * (outputDiagram->GetNumberOfCells() - 1)};
  outputDiagram->InsertNextCell(VTK_LINE, 2, diag.data());

  // add diagonal data
  auto cellData = outputDiagram->GetCellData();
  auto pairIdentifierScalars = vtkIntArray::SafeDownCast(
    cellData->GetArray(ttk::PersistencePairIdentifierName));
  auto extremumIndexScalars = vtkIntArray::SafeDownCast(
    cellData->GetArray(ttk::PersistencePairTypeName));
  auto persistenceScalars = cellData->GetArray(ttk::PersistenceName);
  auto birthScalars = cellData->GetArray(ttk::PersistenceBirthName);
  auto isFinite = cellData->GetArray(ttk::PersistenceIsFinite);

  pairIdentifierScalars->InsertNextTuple1(-1);
  extremumIndexScalars->InsertNextTuple1(-1);
  isFinite->InsertNextTuple1(0);
  // 2 * persistence of min-max pair
  persistenceScalars->InsertNextTuple1(2 * persistenceScalars->GetTuple1(0));
  // birth == death == 0
  birthScalars->InsertNextTuple1(0);

  // don't forget to forward the Field Data
  outputDiagram->GetFieldData()->ShallowCopy(inputDiagram->GetFieldData());

  dbg.printMsg("Projected Persistence Diagram back to 2D", 1.0,
               tm.getElapsedTime(), dbg.getThreadNumber());

  return 0;
}

int TranslateDiagram(vtkUnstructuredGrid *const diagram,
                     const std::array<double, 3> &trans) {

  vtkNew<vtkUnstructuredGrid> tmp{};
  tmp->ShallowCopy(diagram);

  vtkNew<vtkTransform> tr{};
  tr->Translate(trans.data());

  vtkNew<vtkTransformFilter> trf{};
  trf->SetTransform(tr);
  trf->SetInputData(tmp);
  trf->Update();

  diagram->ShallowCopy(trf->GetOutputDataObject(0));

  return 0;
}

int ResetDiagramPosition(vtkUnstructuredGrid *const diagram,
                         const ttk::Debug &dbg) {

  const bool embedded
    = diagram->GetPointData()->GetArray(ttk::PersistenceCoordinatesName)
      == nullptr;

  if(embedded) {
    dbg.printWrn("Cannot reset embedded diagram position");
    return 1;
  }

  // position of first point in diagram
  std::array<double, 3> pos{};
  diagram->GetPoint(diagram->GetCell(0)->GetPointId(0), pos.data());

  // birth value of the first cell in the diagram
  const auto firstBirth{
    diagram->GetCellData()->GetArray(ttk::PersistenceBirthName)->GetTuple1(0)};

  if((pos[0] != pos[1] && pos[0] != firstBirth) || pos[2] != 0) {
    vtkNew<vtkUnstructuredGrid> tmp{};
    tmp->ShallowCopy(diagram);

    vtkNew<vtkTransform> tr{};
    tr->Translate(firstBirth - pos[0], firstBirth - pos[1], -pos[2]);

    vtkNew<vtkTransformFilter> trf{};
    trf->SetTransform(tr);
    trf->SetInputData(tmp);
    trf->Update();

    diagram->ShallowCopy(trf->GetOutputDataObject(0));

    dbg.printMsg("Diagram reset to initial position");
  }

  return 0;
}
