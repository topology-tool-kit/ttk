#include <ttkPersistenceDiagramDistanceMatrix.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkPersistenceDiagramDistanceMatrix);

ttkPersistenceDiagramDistanceMatrix::ttkPersistenceDiagramDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceDiagramDistanceMatrix::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceDiagramDistanceMatrix::updateProgress(const float &progress) {
  {
    std::stringstream msg;
    msg << "[ttkPersistenceDiagramDistanceMatrix] " << progress * 100
        << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDistanceMatrix::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  // Get input data
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  if(blocks != nullptr) {
    inputDiagrams.resize(blocks->GetNumberOfBlocks());
    for(size_t i = 0; i < inputDiagrams.size(); ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    }
  }

  const int numInputs = inputDiagrams.size();

  // Set output
  auto diagramsDistTable = vtkTable::GetData(outputVector);

  std::vector<std::vector<ttk::DiagramTuple>> intermediateDiagrams(numInputs);

  double max_dimension_total = 0.0;
  for(int i = 0; i < numInputs; i++) {
    double max_dimension
      = getPersistenceDiagram(intermediateDiagrams[i], inputDiagrams[i]);
    if(max_dimension_total < max_dimension) {
      max_dimension_total = max_dimension;
    }
  }

  ttk::PersistenceDiagramDistanceMatrix worker{};
  worker.setWrapper(this);

  worker.setWasserstein(WassersteinMetric);
  worker.setPairTypeClustering(PairTypeClustering);
  worker.setAlpha(Alpha);
  worker.setDeltaLim(DeltaLim);
  worker.setLambda(Lambda);
  worker.setUseFullDiagrams(UseFullDiagrams);
  const auto diagramsDistMat = worker.execute(intermediateDiagrams);

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy diagrams distance matrix to output
  for(size_t i = 0; i < diagramsDistMat.size(); ++i) {
    std::string name{"Diagram"};
    zeroPad(name, diagramsDistMat.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numInputs);
    col->SetName(name.c_str());
    for(size_t j = 0; j < diagramsDistMat[i].size(); ++j) {
      col->SetTuple1(j, diagramsDistMat[i][j]);
    }
    diagramsDistTable->AddColumn(col);
  }

  // aggregate input field data
  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputDiagrams[0]->GetFieldData());
  fd->SetNumberOfTuples(inputDiagrams.size());
  for(size_t i = 0; i < inputDiagrams.size(); ++i) {
    fd->SetTuple(i, 0, inputDiagrams[i]->GetFieldData());
  }

  // copy input field data to output row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    diagramsDistTable->AddColumn(fd->GetAbstractArray(i));
  }

  {
    std::stringstream msg;
    msg << "[ttkPersistenceDiagramDistanceMatrix] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}

double ttkPersistenceDiagramDistanceMatrix::getPersistenceDiagram(
  std::vector<ttk::DiagramTuple> &diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_) {

  const auto pd = CTPersistenceDiagram_->GetPointData();
  const auto cd = CTPersistenceDiagram_->GetCellData();
  const auto points = CTPersistenceDiagram_->GetPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(pd == nullptr || points == nullptr) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto vertexIdentifierScalars
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pd->GetArray("CriticalType"));
  const auto pairIdentifierScalars
    = vtkIntArray::SafeDownCast(cd->GetArray("PairIdentifier"));
  const auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cd->GetArray("PairType"));
  const auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cd->GetArray("Persistence"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Death"));
  const auto critCoordinates
    = vtkFloatArray::SafeDownCast(pd->GetArray("Coordinates"));

  const bool embed = birthScalars != nullptr && deathScalars != nullptr;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!embed && critCoordinates == nullptr) {
    // missing data
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
      pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  diagram.resize(pairingsSize + 1);
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize - 1; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    std::array<double, 3> coordsBirth{}, coordsDeath{};

    const auto i0 = 2 * i;
    const auto i1 = 2 * i + 1;

    double birth, death;

    if(embed) {
      points->GetPoint(i0, coordsBirth.data());
      points->GetPoint(i1, coordsDeath.data());
      birth = birthScalars->GetValue(i0);
      death = deathScalars->GetValue(i1);
    } else {
      critCoordinates->GetTuple(i0, coordsBirth.data());
      critCoordinates->GetTuple(i1, coordsDeath.data());
      birth = points->GetPoint(i0)[0];
      death = points->GetPoint(i1)[1];
    }

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
      if(pairIdentifier == 0) {
        max_dimension = persistence;

        diagram[0] = std::make_tuple(
          vertexId1, ttk::CriticalType::Local_minimum, vertexId2,
          ttk::CriticalType::Saddle1, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);
        diagram[pairingsSize] = std::make_tuple(
          vertexId1, ttk::CriticalType::Saddle1, vertexId2,
          ttk::CriticalType::Local_maximum, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);

      } else {
        diagram[pairIdentifier] = std::make_tuple(
          vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
          persistence, pairType, birth, coordsBirth[0], coordsBirth[1],
          coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
      }
    }
    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKPersistenceDiagramDistanceMatrix] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramDistanceMatrix] Missed " << nbNonCompact
          << " pairs due to non-compactness." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  return max_dimension;
}
