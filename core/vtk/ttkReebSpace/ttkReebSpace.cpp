#include <vtkNew.h>

#include <ttkMacros.h>
#include <ttkReebSpace.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkReebSpace);
ttkReebSpace::ttkReebSpace() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(4);
}

int ttkReebSpace::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkReebSpace::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) { // 0-sheets, corners of jacobi set segments
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 1) { // 1-sheets, jacobi sets
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) { // 2-sheets, fiber surfaces of jacobi sets
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <class dataTypeU, class dataTypeV>
int ttkReebSpace::baseCall(vtkDataSet *input,
                           vtkDataArray *uField,
                           vtkDataArray *offsetFieldU,
                           vtkDataArray *vField,
                           vtkDataArray *offsetFieldV) {

  bool VaryingValues
    = (uField->GetMTime() > GetMTime())
      || ((offsetFieldU) && (offsetFieldU->GetMTime() > GetMTime()))
      || (vField->GetMTime() > GetMTime())
      || ((offsetFieldV) && (offsetFieldV->GetMTime() > GetMTime()))
      || (this->setRangeDrivenOctree(UseOctreeAcceleration));

  // first time or the values changed
  this->setSosOffsetsU(&sosOffsetsU_);
  this->setSosOffsetsV(&sosOffsetsV_);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(!triangulation)
    return -1;

  bool VaryingTriangulation = false;
  if(triangulation->isEmpty())
    VaryingTriangulation = true;

  this->preconditionTriangulation<dataTypeU, dataTypeV>(triangulation);

  this->printMsg("U-component: `" + std::string{uField->GetName()} + "'");
  this->printMsg("V-component: `" + std::string{vField->GetName()} + "'");

  // go!
  if((this->empty()) || (VaryingValues) || (VaryingTriangulation)) {
    this->printMsg("Starting computation");
    this->execute(static_cast<dataTypeU *>(ttkUtils::GetVoidPointer(uField)),
                  static_cast<dataTypeV *>(ttkUtils::GetVoidPointer(vField)),
                  *triangulation->getData());
  }

  if(SimplificationThreshold > 0) {
    this->simplify(static_cast<dataTypeU *>(ttkUtils::GetVoidPointer(uField)),
                   static_cast<dataTypeV *>(ttkUtils::GetVoidPointer(vField)),
                   *triangulation->getData(), SimplificationThreshold,
                   (ReebSpace::SimplificationCriterion)SimplificationCriterion);
  }

  Modified();

  return 0;
}
int ttkReebSpace::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto sheet0 = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto sheet1 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto sheet2 = vtkUnstructuredGrid::GetData(outputVector, 2);
  auto sheet3 = vtkDataSet::GetData(outputVector, 3);

  // check data components
  if(Ucomponent.length()) {
    std::string oldName = Ucomponent;
    if(uComponent_) {
      oldName = uComponent_->GetName();
    }
    uComponent_ = input->GetPointData()->GetArray(Ucomponent.data());
    if(oldName != Ucomponent)
      uComponent_->Modified();
  } else {
    // default
    uComponent_ = input->GetPointData()->GetArray(UcomponentId);
  }
  if(!uComponent_)
    return -1;

  if(Vcomponent.length()) {
    std::string oldName = Vcomponent;
    if(vComponent_) {
      oldName = vComponent_->GetName();
    }
    vComponent_ = input->GetPointData()->GetArray(Vcomponent.data());
    if(oldName != Vcomponent)
      vComponent_->Modified();
  } else {
    // default
    vComponent_ = input->GetPointData()->GetArray(VcomponentId);
  }
  if(!vComponent_)
    return -2;

  if(ForceInputOffsetScalarField) {
    if(OffsetFieldU.length()) {

      std::string oldName = OffsetFieldU;
      if(offsetFieldU_) {
        oldName = offsetFieldU_->GetName();
      }
      offsetFieldU_ = input->GetPointData()->GetArray(OffsetFieldU.data());
      if(oldName != OffsetFieldU)
        offsetFieldU_->Modified();

      if(offsetFieldU_) {
        sosOffsetsU_.resize(offsetFieldU_->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldU_->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU_->GetTuple1(i);
        }
      }
    }
    if(OffsetFieldV.length()) {

      std::string oldName = OffsetFieldV;
      if(offsetFieldV_) {
        oldName = offsetFieldV_->GetName();
      }
      offsetFieldV_ = input->GetPointData()->GetArray(OffsetFieldV.data());
      if(oldName != OffsetFieldV)
        offsetFieldV_->Modified();

      if(offsetFieldV_) {
        sosOffsetsV_.resize(offsetFieldV_->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldV_->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV_->GetTuple1(i);
        }
      }
    }
  } else {
    if(input->GetPointData()->GetArray(ttk::OffsetFieldUName)) {
      offsetFieldU_ = input->GetPointData()->GetArray(ttk::OffsetFieldUName);
      if(offsetFieldU_) {
        sosOffsetsU_.resize(offsetFieldU_->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldU_->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU_->GetTuple1(i);
        }
      }
    }
    if(input->GetPointData()->GetArray(ttk::OffsetFieldVName)) {
      offsetFieldV_ = input->GetPointData()->GetArray(ttk::OffsetFieldVName);
      if(offsetFieldV_) {
        sosOffsetsV_.resize(offsetFieldV_->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldV_->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV_->GetTuple1(i);
        }
      }
    }
  }

  // set the Reeb space functor
  switch(vtkTemplate2PackMacro(
    uComponent_->GetDataType(), vComponent_->GetDataType())) {
    vtkTemplate2Macro((baseCall<VTK_T1, VTK_T2>(
      input, uComponent_, offsetFieldU_, vComponent_, offsetFieldV_)));
  }

  // prepare the output
  this->printMsg("Preparing the VTK output...");

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(!triangulation)
    return -3;

  // 0-sheets -
  // Optional additional fields:
  // PointData; u, v, vertexId, type, sheetId
  const auto sheet0segmentation = this->get0sheetSegmentation();
  int vertexNumber = 0;
  for(size_t i = 0; i < sheet0segmentation->size(); i++) {
    int sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1) {
      const ReebSpace::Sheet0 *sheet = this->get0sheet(sheet0Id);
      if(sheet != nullptr && !sheet->pruned_) {
        vertexNumber++;
      }
    }
  }
  vtkNew<vtkPoints> sheet0Points{};
  if(!sheet0->GetPoints())
    sheet0->SetPoints(sheet0Points);

  sheet0->GetPoints()->SetNumberOfPoints(vertexNumber);

  vtkNew<vtkDoubleArray> vertexScalarsU{};
  vtkNew<vtkDoubleArray> vertexScalarsV{};
  if(ZeroSheetValue) {
    vertexScalarsU->SetNumberOfTuples(vertexNumber);
    vertexScalarsU->SetName(uComponent_->GetName());
    vertexScalarsV->SetNumberOfTuples(vertexNumber);
    vertexScalarsV->SetName(vComponent_->GetName());
  } else {
    sheet0->GetPointData()->RemoveArray(uComponent_->GetName());
    sheet0->GetPointData()->RemoveArray(vComponent_->GetName());
  }

  vtkNew<ttkSimplexIdTypeArray> vertexIds{};
  if(ZeroSheetVertexId) {
    vertexIds->SetNumberOfTuples(vertexNumber);
    vertexIds->SetName(ttk::VertexScalarFieldName);
  } else {
    sheet0->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  vtkNew<vtkCharArray> vertexTypes{};
  if(ZeroSheetType) {
    vertexTypes->SetNumberOfTuples(vertexNumber);
    vertexTypes->SetName("SheetType");
  } else {
    sheet0->GetPointData()->RemoveArray("SheetType");
  }

  vtkNew<ttkSimplexIdTypeArray> vertexSheetId{};
  if(ZeroSheetId) {
    vertexSheetId->SetNumberOfTuples(vertexNumber);
    vertexSheetId->SetName("0-SheetId");
  } else {
    sheet0->GetPointData()->RemoveArray("0-SheetId");
  }

  vertexNumber = 0;
  double *p = nullptr;
  for(size_t i = 0; i < sheet0segmentation->size(); i++) {
    int sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1) {

      const ReebSpace::Sheet0 *sheet = this->get0sheet(sheet0Id);

      if(sheet != nullptr && !sheet->pruned_) {
        p = input->GetPoint(i);

        sheet0->GetPoints()->SetPoint(vertexNumber, p);

        if(ZeroSheetId) {
          vertexSheetId->SetTuple1(vertexNumber, (*sheet0segmentation)[i]);
        }
        if(ZeroSheetVertexId) {
          vertexIds->SetTuple1(vertexNumber, i);
        }
        if(ZeroSheetValue) {
          double u, v;

          uComponent_->GetTuple(i, &u);
          vComponent_->GetTuple(i, &v);

          vertexScalarsU->SetTuple1(vertexNumber, u);
          vertexScalarsV->SetTuple1(vertexNumber, v);
        }
        if(ZeroSheetType) {
          const ReebSpace::Sheet0 *sht0
            = this->get0sheet((*sheet0segmentation)[i]);
          vertexTypes->SetTuple1(vertexNumber, sht0->type_);
        }

        vertexNumber++;
      }
    }
  }
  if(ZeroSheetId)
    sheet0->GetPointData()->AddArray(vertexSheetId);
  if(ZeroSheetVertexId)
    sheet0->GetPointData()->AddArray(vertexIds);
  if(ZeroSheetValue) {
    sheet0->GetPointData()->AddArray(vertexScalarsU);
    sheet0->GetPointData()->AddArray(vertexScalarsV);
  }
  if(ZeroSheetType)
    sheet0->GetPointData()->AddArray(vertexTypes);

  // 1-sheets
  // Optional additional fields:
  // PointData: u, v, vertexId,
  // CellData: edgeId, type, sheetId
  const auto sheet1segmentation = this->get1sheetSegmentation();

  vtkNew<vtkPoints> sheet1Points{};
  sheet1->SetPoints(sheet1Points);

  vtkNew<vtkCellArray> sheet1Edges{};
  vtkNew<vtkDoubleArray> edgeScalarsU{};
  vtkNew<vtkDoubleArray> edgeScalarsV{};

  if(OneSheetValue) {
    edgeScalarsU->SetName(uComponent_->GetName());
    edgeScalarsV->SetName(vComponent_->GetName());
  } else {
    sheet1->GetPointData()->RemoveArray(uComponent_->GetName());
    sheet1->GetPointData()->RemoveArray(vComponent_->GetName());
  }

  vtkNew<ttkSimplexIdTypeArray> edgeVertexIds{};
  if(OneSheetVertexId) {
    edgeVertexIds->SetName(ttk::VertexScalarFieldName);
  } else {
    sheet1->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  vtkNew<vtkIntArray> edgeType{};
  if(OneSheetType) {
    edgeType->SetName("EdgeType");
  } else {
    sheet1->GetCellData()->RemoveArray("EdgeType");
  }

  vtkNew<ttkSimplexIdTypeArray> edgeIds{};
  if(OneSheetEdgeId) {
    edgeIds->SetName("EdgeIds");
  } else {
    sheet1->GetCellData()->RemoveArray("EdgeIds");
  }

  vtkNew<ttkSimplexIdTypeArray> edgeSheetIds{};
  if(OneSheetId) {
    edgeSheetIds->SetName("1-SheetId");
  } else {
    sheet1->GetCellData()->RemoveArray("1-SheetId");
  }

  vertexNumber = 0;
  double p0[3], p1[3];
  vtkNew<vtkIdList> idList{};
  idList->SetNumberOfIds(2);
  const auto edgeTypes = this->getEdgeTypes();

  for(size_t i = 0; i < sheet1segmentation->size(); i++) {

    int sheet1Id = (*sheet1segmentation)[i];

    if(sheet1Id != -1) {

      const ReebSpace::Sheet1 *sheet = this->get1sheet(sheet1Id);

      if((sheet) && (!sheet->pruned_)) {

        int vertexId0 = -1, vertexId1 = -1;
        triangulation->getEdgeVertex(i, 0, vertexId0);
        triangulation->getEdgeVertex(i, 1, vertexId1);

        input->GetPoint(vertexId0, p0);
        input->GetPoint(vertexId1, p1);

        sheet1->GetPoints()->InsertNextPoint(p0);
        sheet1->GetPoints()->InsertNextPoint(p1);

        if(OneSheetValue) {
          double u, v;
          uComponent_->GetTuple(vertexId0, &u);
          vComponent_->GetTuple(vertexId0, &v);

          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);

          uComponent_->GetTuple(vertexId1, &u);
          vComponent_->GetTuple(vertexId1, &v);

          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);
        }

        if(OneSheetVertexId) {
          edgeVertexIds->InsertNextTuple1(vertexId0);
          edgeVertexIds->InsertNextTuple1(vertexId1);
        }

        idList->SetId(0, vertexNumber);
        idList->SetId(1, vertexNumber + 1);
        sheet1Edges->InsertNextCell(idList);
        vertexNumber += 2;

        if(OneSheetEdgeId) {
          edgeIds->InsertNextTuple1(i);
        }
        if(OneSheetType) {
          edgeType->InsertNextTuple1((*edgeTypes)[i]);
        }
        if(OneSheetId) {
          edgeSheetIds->InsertNextTuple1((*sheet1segmentation)[i]);
        }
      }
    }
  }

  sheet1->SetCells(VTK_LINE, sheet1Edges);

  if(OneSheetValue) {
    sheet1->GetPointData()->AddArray(edgeScalarsU);
    sheet1->GetPointData()->AddArray(edgeScalarsV);
  }
  if(OneSheetVertexId) {
    sheet1->GetPointData()->AddArray(edgeVertexIds);
  }
  if(OneSheetId) {
    sheet1->GetCellData()->AddArray(edgeSheetIds);
  }
  if(OneSheetEdgeId) {
    sheet1->GetCellData()->AddArray(edgeIds);
  }
  if(OneSheetType) {
    sheet1->GetCellData()->AddArray(edgeType);
  }

  // 2-sheets
  // optional fields:
  // pointdata: twoSheetValues, twoSheetParameterization
  if(TwoSheets) {
    const std::vector<ttk::FiberSurface::Vertex> *vertexList
      = this->getFiberSurfaceVertices();

    vtkNew<vtkPoints> sheet2Points{};
    sheet2Points->SetNumberOfPoints(vertexList->size());
    sheet2->SetPoints(sheet2Points);

    vtkNew<vtkDoubleArray> triangleScalarsU{};
    vtkNew<vtkDoubleArray> triangleScalarsV{};

    if(TwoSheetValue) {
      triangleScalarsU->SetName(uComponent_->GetName());
      triangleScalarsU->SetNumberOfTuples(vertexList->size());
      triangleScalarsV->SetName(vComponent_->GetName());
      triangleScalarsV->SetNumberOfTuples(vertexList->size());
    } else {
      sheet2->GetPointData()->RemoveArray(uComponent_->GetName());
      sheet2->GetPointData()->RemoveArray(vComponent_->GetName());
    }

    vtkNew<vtkDoubleArray> triangleParameterization{};
    if(TwoSheetParameterization) {
      triangleParameterization->SetName("EdgeParameterization");
      triangleParameterization->SetNumberOfTuples(vertexList->size());
    } else {
      sheet2->GetPointData()->RemoveArray("EdgeParameterization");
    }

    int sheet2TriangleNumber = 0;
    for(int i = 0; i < this->getNumberOf2sheets(); i++) {
      const ReebSpace::Sheet2 *sheet = this->get2sheet(i);

      if(sheet != nullptr && !sheet->pruned_) {
        for(size_t j = 0; j < sheet->triangleList_.size(); j++) {
          sheet2TriangleNumber += sheet->triangleList_[j].size();
        }
      }
    }

    vtkNew<vtkCellArray> sheet2Triangles{};

    // celldata: twoSheetId, twoSheetEdgeId, twoSheetTetId
    vtkNew<ttkSimplexIdTypeArray> triangleSheetIds{};
    if(TwoSheetId) {
      triangleSheetIds->SetName("2-SheetId");
      triangleSheetIds->SetNumberOfTuples(sheet2TriangleNumber);
    } else {
      sheet2->GetCellData()->RemoveArray("2-SheetId");
    }

    vtkNew<ttkSimplexIdTypeArray> triangleEdgeIds{};
    if(TwoSheetEdgeId) {
      triangleEdgeIds->SetName("EdgeIds");
      triangleEdgeIds->SetNumberOfTuples(sheet2TriangleNumber);
    } else {
      sheet2->GetCellData()->RemoveArray("EdgeIds");
    }

    vtkNew<vtkIntArray> triangleEdgeType{};
    if(TwoSheetEdgeType) {
      triangleEdgeType->SetName("EdgeType");
      triangleEdgeType->SetNumberOfTuples(sheet2TriangleNumber);
    } else {
      sheet2->GetCellData()->RemoveArray("EdgeType");
    }

    vtkNew<ttkSimplexIdTypeArray> triangleTetIds{};
    if(TwoSheetTetId) {
      triangleTetIds->SetName("TetIds");
      triangleTetIds->SetNumberOfTuples(sheet2TriangleNumber);
    } else {
      sheet2->GetCellData()->RemoveArray("TetIds");
    }

    vtkNew<ttkSimplexIdTypeArray> triangleCaseIds{};
    if(TwoSheetCaseId) {
      triangleCaseIds->SetName("CaseIds");
      triangleCaseIds->SetNumberOfTuples(sheet2TriangleNumber);
    } else {
      sheet2->GetCellData()->RemoveArray("CaseIds");
    }

    for(size_t i = 0; i < vertexList->size(); i++) {
      sheet2->GetPoints()->SetPoint(i, (*vertexList)[i].p_[0],
                                    (*vertexList)[i].p_[1],
                                    (*vertexList)[i].p_[2]);

      if(TwoSheetValue) {
        triangleScalarsU->SetTuple1(i, (*vertexList)[i].uv_.first);
        triangleScalarsV->SetTuple1(i, (*vertexList)[i].uv_.second);
      }
      if(TwoSheetParameterization) {
        triangleParameterization->SetTuple1(i, (*vertexList)[i].t_);
      }
    }
    if(TwoSheetValue) {
      sheet2->GetPointData()->AddArray(triangleScalarsU);
      sheet2->GetPointData()->AddArray(triangleScalarsV);
    }
    if(TwoSheetParameterization) {
      sheet2->GetPointData()->AddArray(triangleParameterization);
    }

    int triangleNumber = 0;
    idList->SetNumberOfIds(3);
    for(int i = 0; i < this->getNumberOf2sheets(); i++) {
      const ReebSpace::Sheet2 *sht2 = this->get2sheet(i);

      if(sht2 != nullptr && !sht2->pruned_) {
        for(size_t j = 0; j < sht2->triangleList_.size(); j++) {

          for(size_t k = 0; k < sht2->triangleList_[j].size(); k++) {

            for(int l = 0; l < 3; l++) {
              idList->SetId(l, sht2->triangleList_[j][k].vertexIds_[l]);
            }

            sheet2Triangles->InsertNextCell(idList);

            if(TwoSheetId) {
              triangleSheetIds->SetTuple1(triangleNumber, i);
            }

            if(TwoSheetEdgeId) {
              const ReebSpace::Sheet1 *sht1 = this->get1sheet(sht2->sheet1Id_);
              triangleEdgeIds->SetTuple1(triangleNumber, sht1->edgeList_[j]);
            }

            if(TwoSheetEdgeType) {
              auto polygonEdgeId = sht2->triangleList_[j][k].polygonEdgeId_;
              auto edgeId = this->getJacobi2Edge(polygonEdgeId);
              triangleEdgeType->SetTuple1(triangleNumber, (*edgeTypes)[edgeId]);
            }

            if(TwoSheetTetId) {
              triangleTetIds->SetTuple1(
                triangleNumber, sht2->triangleList_[j][k].tetId_);
            }
            if(TwoSheetCaseId) {
              triangleCaseIds->SetTuple1(
                triangleNumber, sht2->triangleList_[j][k].caseId_);
            }

            triangleNumber++;
          }
        }
      }
    }
    sheet2->SetCells(VTK_TRIANGLE, sheet2Triangles);

    if(TwoSheetId) {
      sheet2->GetCellData()->AddArray(triangleSheetIds);
    }
    if(TwoSheetEdgeId) {
      sheet2->GetCellData()->AddArray(triangleEdgeIds);
    }
    if(TwoSheetEdgeType) {
      sheet2->GetCellData()->AddArray(triangleEdgeType);
    }
    if(TwoSheetTetId) {
      sheet2->GetCellData()->AddArray(triangleTetIds);
    }
    if(TwoSheetCaseId) {
      sheet2->GetCellData()->AddArray(triangleCaseIds);
    }
  }

  // now take care of the 3 sheets
  //   vector<float> *triangulationPoints
  //     = reebSpace_.getSheetTriangulationPoints();
  //   vector<long long int> *triangulationCells
  //     = reebSpace_.getSheetTriangulationCells();
  //
  //   vtkSmartPointer<vtkPoints> sheet3Points =
  //   vtkSmartPointer<vtkPoints>::New(); vtkSmartPointer<vtkFloatArray>
  //   pointData =
  //     vtkSmartPointer<vtkFloatArray>::New();
  //   pointData->SetNumberOfComponents(3);
  //   ttkUtils::SetVoidArray(
  //     pointData, triangulationPoints->data(), triangulationPoints->size(),
  //     1);
  //   sheet3Points->SetData(pointData);
  //   sheet3->SetPoints(sheet3Points);
  //
  //   vtkSmartPointer<vtkCellArray> sheet3Cells
  //     = vtkSmartPointer<vtkCellArray>::New();
  //   vtkSmartPointer<ttkSimplexIdTypeArray> idArray
  //     = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  //   ttkUtils::SetVoidArray(
  //     idArray, triangulationCells->data(), triangulationCells->size(), 1);
  //   sheet3Cells->SetCells(triangulationCells->size()/5, idArray);
  //   sheet3->SetCells(VTK_TETRA, sheet3Cells);

  // now take care of the 3 sheets
  sheet3->ShallowCopy(input);
  const auto vertex3sheets = this->get3sheetVertexSegmentation();

  vtkNew<ttkSimplexIdTypeArray> vertexNumberField{};
  vtkNew<ttkSimplexIdTypeArray> tetNumberField{};

  if(ThreeSheetTetNumber) {
    tetNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    tetNumberField->SetName("3-SheetTetNumber");
    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      const ReebSpace::Sheet3 *sht3 = this->get3sheet((*vertex3sheets)[i]);
      if((sht3) && (!sht3->pruned_))
        tetNumberField->SetTuple1(i, sht3->tetList_.size());
      else
        tetNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(tetNumberField);
  } else {
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }

  if(ThreeSheetVertexNumber) {
    vertexNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    vertexNumberField->SetName("3-SheetVertexNumber");
    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      const ReebSpace::Sheet3 *sht3 = this->get3sheet((*vertex3sheets)[i]);
      if((sht3) && (!sht3->pruned_))
        vertexNumberField->SetTuple1(i, sht3->vertexList_.size());
      else
        vertexNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(vertexNumberField);
  } else {
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }

  vtkNew<vtkDoubleArray> domainVolume{};
  if(ThreeSheetDomainVolume) {
    domainVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    domainVolume->SetName("3-SheetDomainVolume");

    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      const ReebSpace::Sheet3 *sht3 = this->get3sheet((*vertex3sheets)[i]);
      if((sht3) && (!sht3->pruned_)) {
        domainVolume->SetTuple1(i, sht3->domainVolume_);
      } else {
        domainVolume->SetTuple1(i, 0);
      }
    }

    sheet3->GetPointData()->AddArray(domainVolume);
  } else {
    sheet3->GetPointData()->RemoveArray("3-SheetDomainVolume");
  }

  vtkNew<vtkDoubleArray> rangeArea{};
  if(ThreeSheetDomainVolume) {
    rangeArea->SetNumberOfTuples(input->GetNumberOfPoints());
    rangeArea->SetName("3-SheetRangeArea");

    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      const ReebSpace::Sheet3 *sht3 = this->get3sheet((*vertex3sheets)[i]);
      if((sht3) && (!sht3->pruned_)) {
        rangeArea->SetTuple1(i, sht3->rangeArea_);
      } else {
        rangeArea->SetTuple1(i, 0);
      }
    }

    sheet3->GetPointData()->AddArray(rangeArea);
  } else {
    sheet3->GetPointData()->RemoveArray("3-SheetRangeArea");
  }

  vtkNew<vtkDoubleArray> hyperVolume{};
  if(ThreeSheetDomainVolume) {
    hyperVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    hyperVolume->SetName("3-SheetHyperVolume");

    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      const ReebSpace::Sheet3 *sht3 = this->get3sheet((*vertex3sheets)[i]);
      if((sht3) && (!sht3->pruned_)) {
        hyperVolume->SetTuple1(i, sht3->hyperVolume_);
      } else {
        hyperVolume->SetTuple1(i, 0);
      }
    }

    sheet3->GetPointData()->AddArray(hyperVolume);
  } else {
    sheet3->GetPointData()->RemoveArray("3-SheetHyperVolume");
  }

  vtkNew<ttkSimplexIdTypeArray> vertexSegmentation{};
  vertexSegmentation->SetName("3-SheetId");
  vertexSegmentation->SetNumberOfTuples(input->GetNumberOfPoints());
  for(int i = 0; i < input->GetNumberOfPoints(); i++) {
    const ReebSpace::Sheet3 *sheet = this->get3sheet((*vertex3sheets)[i]);
    if(sheet) {
      vertexSegmentation->SetTuple1(i, sheet->simplificationId_);
    } else {
      vertexSegmentation->SetTuple1(i, (*vertex3sheets)[i]);
    }
  }
  sheet3->GetPointData()->AddArray(vertexSegmentation);

  const auto tet3sheets = this->get3sheetTetSegmentation();
  vtkNew<ttkSimplexIdTypeArray> tetSegmentation{};
  tetSegmentation->SetName("3-SheetId");
  tetSegmentation->SetNumberOfTuples(input->GetNumberOfCells());
  for(int i = 0; i < input->GetNumberOfCells(); i++) {
    const ReebSpace::Sheet3 *sht3 = this->get3sheet((*tet3sheets)[i]);
    if(sht3) {
      tetSegmentation->SetTuple1(i, sht3->simplificationId_);
    } else {
      tetSegmentation->SetTuple1(i, (*tet3sheets)[i]);
    }
  }
  sheet3->GetCellData()->AddArray(tetSegmentation);

  return 0;
}
