#include <ttkPersistenceDiagramDictionaryDecoding.h>

#include <vtkAbstractArray.h>
#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkVariantArray.h>

#include <ttkPersistenceDiagramUtils.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagramDictionaryDecoding);

ttkPersistenceDiagramDictionaryDecoding::
  ttkPersistenceDiagramDictionaryDecoding() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

int ttkPersistenceDiagramDictionaryDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

int ttkPersistenceDiagramDictionaryDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

int ttkPersistenceDiagramDictionaryDecoding::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  const auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0]);
  const auto weightsVTK = vtkTable::GetData(inputVector[1]);

  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  if(blocks != nullptr) {
    int numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    }
  }

  const size_t nDiags = inputDiagrams.size();

  std::vector<ttk::DiagramType> dictDiagrams(nDiags);
  for(size_t i = 0; i < nDiags; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu;
    vtu->DeepCopy(inputDiagrams[i]);
    auto &atom = dictDiagrams[i];
    const auto ret = VTUToDiagram(atom, vtu, *this);
    if(ret != 0) {
      this->printWrn("Could not read Persistence Diagram");
    }
  }

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  vtkNew<vtkTable> temp;
  temp->DeepCopy(weightsVTK);
  std::vector<vtkDataArray *> inputWeights;
  int numWeights = temp->GetNumberOfRows();
  for(int i = 0; i < temp->GetNumberOfColumns(); ++i) {
    std::cout << temp->GetColumnName(i) << "\n";
  }

  if(weightsVTK != nullptr) {
    inputWeights.resize(nDiags);
    for(size_t i = 0; i < nDiags; ++i) {
      std::string name{"Atom"};
      zeroPad(name, nDiags, i);
      inputWeights[i]
        = vtkDataArray::SafeDownCast(temp->GetColumnByName(name.c_str()));
    }
  }

  std::vector<std::vector<double>> vectorWeights(numWeights);
  for(int i = 0; i < numWeights; ++i) {
    std::vector<double> &t1 = vectorWeights[i];
    for(size_t j = 0; j < nDiags; ++j) {
      double weight = inputWeights[j]->GetTuple1(i);
      t1.push_back(weight);
    }
  }
  std::vector<ttk::DiagramType> Barycenters(numWeights);

  if(!ComputePoints) {
    this->execute(dictDiagrams, vectorWeights, Barycenters);
  }

  auto outputDgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto outputCoordinates = vtkTable::GetData(outputVector, 1);
  outputDgm->SetNumberOfBlocks(numWeights);
  outputCoordinates->SetNumberOfRows(numWeights);

  outputDiagrams(outputDgm, outputCoordinates, Barycenters, dictDiagrams,
                 weightsVTK, vectorWeights, Spacing, 1);

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation

  // make a SHALLOW copy of the input
  // outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  // outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}

void ttkPersistenceDiagramDictionaryDecoding::outputDiagrams(
  vtkMultiBlockDataSet *output,
  vtkTable *outputCoordinates,
  const std::vector<ttk::DiagramType> &diags,
  std::vector<ttk::DiagramType> &atoms,
  vtkTable *weightsVTK,
  const std::vector<std::vector<double>> &weights,
  const double spacing,
  const double maxPersistence) const {

  const auto nDiags = diags.size();
  const auto nAtoms = atoms.size();

  ttk::SimplexId n_existing_blocks = ShowAtoms ? nAtoms : 0;

  output->SetNumberOfBlocks(nDiags + n_existing_blocks);
  std::vector<std::array<double, 3>> coords(nAtoms);
  std::vector<std::array<double, 3>> trueCoords(nAtoms);
  std::vector<double> xVector(nDiags);
  std::vector<double> yVector(nDiags);
  std::vector<double> zVector(nDiags, 0.);
  vtkNew<vtkDoubleArray> dummy{};

  computeAtomsCoordinates(atoms, weights, coords, trueCoords, xVector, yVector,
                          zVector, spacing, nAtoms);

  if(nAtoms == 2) {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i][0];
        double Y = coords[i][1];
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else if(nAtoms == 3) {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i][0];
        double Y = coords[i][1];
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        const auto angle
          = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nAtoms);
        double X = spacing * maxPersistence * std::cos(angle);
        double Y = spacing * maxPersistence * std::sin(angle);
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }
  }

  int numDiags = diags.size();

  for(int i = 0; i < numDiags; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    if(!ComputePoints) {
      DiagramToVTU(vtu, diags[i], dummy, *this, 3, false);
    }

    double X = 0;
    double Y = 0;
    for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
      X += weights[i][iAtom] * coords[iAtom][0];
      Y += weights[i][iAtom] * coords[iAtom][1];
    }

    TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
    output->SetBlock(i + n_existing_blocks, vtu);
  }

  for(size_t i = 0; i < 3; ++i) {
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(nDiags);
    std::string name;

    if(i == 0) {
      name = "X";
    } else if(i == 1) {
      name = "Y";
    } else {
      name = "Z";
    }
    col->SetName(name.c_str());
    for(size_t j = 0; j < nDiags; ++j) {
      if(i == 0) {
        col->SetValue(j, xVector[j]);
      } else if(i == 1) {
        col->SetValue(j, yVector[j]);
      } else {
        col->SetValue(j, zVector[j]);
      }
    }
    col->Modified();
    outputCoordinates->AddColumn(col);
  }

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  vtkNew<vtkTable> temp;
  temp->DeepCopy(weightsVTK);
  for(int i = 0; i < temp->GetNumberOfColumns(); ++i) {
    int test = 0;
    const auto array = temp->GetColumn(i);

    for(size_t j = 0; j < nDiags; ++j) {
      std::string name{"Atom"};
      zeroPad(name, nDiags, j);
      if(strcmp(name.c_str(), temp->GetColumnName(i)) == 0) {
        test += 1;
      }
    }
    if(test > 0) {
      continue;
    }
    outputCoordinates->AddColumn(array);
  }

  for(size_t i = 0; i < nAtoms; ++i) {
    vtkNew<vtkVariantArray> row{};
    row->SetNumberOfValues(outputCoordinates->GetNumberOfColumns());
    for(int j = 0; j < outputCoordinates->GetNumberOfColumns(); ++j) {
      if(strcmp(outputCoordinates->GetColumnName(j), "X") == 0) {
        row->SetValue(j, trueCoords[i][0]);
      } else if(strcmp(outputCoordinates->GetColumnName(j), "Y") == 0) {
        row->SetValue(j, trueCoords[i][1]);
      } else if(strcmp(outputCoordinates->GetColumnName(j), "Z") == 0) {
        row->SetValue(j, trueCoords[i][2]);
      } else if(strcmp(outputCoordinates->GetColumnName(j), "ClusterID") == 0) {
        row->SetValue(j, -1);
      } else {
        continue;
      }
    }
    row->Modified();
    outputCoordinates->InsertNextRow(row);
  }
}

double ttkPersistenceDiagramDictionaryDecoding::getMaxPersistence(
  const ttk::DiagramType &diagram) const {

  double maxPersistence{0};
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &t = diagram[i];
    const double &pers = t.persistence();
    maxPersistence = std::max(pers, maxPersistence);
  }
  return maxPersistence;
}
