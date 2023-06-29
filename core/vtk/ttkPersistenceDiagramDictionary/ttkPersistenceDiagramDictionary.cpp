#include <ttkMacros.h>
#include <ttkPersistenceDiagramDictionary.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkAlgorithm.h>
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

vtkStandardNewMacro(ttkPersistenceDiagramDictionary);

ttkPersistenceDiagramDictionary::ttkPersistenceDiagramDictionary() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(2);
}

int ttkPersistenceDiagramDictionary::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDictionary::FillOutputPortInformation(
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

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDictionary::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  // Get input data
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto atomBlocks = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // Flat storage for diagrams extracted from blocks
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  std::vector<vtkUnstructuredGrid *> inputAtoms;

  int numInputAtoms = 0;
  if(atomBlocks != nullptr) {
    numInputAtoms = atomBlocks->GetNumberOfBlocks();
    inputAtoms.resize(numInputAtoms);
    for(int i = 0; i < numInputAtoms; ++i) {
      inputAtoms[i]
        = vtkUnstructuredGrid::SafeDownCast(atomBlocks->GetBlock(i));
    }
  }
  if(BackEnd == BACKEND::INPUT_ATOMS) {
    AtomNumber_ = numInputAtoms;
  }
  if(blocks != nullptr) {
    int numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    }
  }

  const int numAtom = this->GetAtomNumber_();
  printMsg("Number of atoms: " + ttk::debug::output::YELLOW
           + ttk::debug::output::UNDERLINED + std::to_string(numAtom)
           + ttk::debug::output::ENDCOLOR + ttk::debug::output::ENDCOLOR);

  // total number of diagrams
  const int nDiags = inputDiagrams.size();

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  // Set output
  auto outputDgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto outputWeights = vtkTable::GetData(outputVector, 1);

  outputDgm->SetNumberOfBlocks(numAtom);

  if(BackEnd == BACKEND::INPUT_ATOMS) {
    for(int i = 0; i < numAtom; ++i) {
      vtkNew<vtkUnstructuredGrid> vtu;
      vtu->DeepCopy(inputAtoms[i]);
      outputDgm->SetBlock(i, vtu);
    }
  } else {
    for(int i = 0; i < numAtom; ++i) {
      vtkNew<vtkUnstructuredGrid> vtu;
      vtu->DeepCopy(inputDiagrams[i]);
      outputDgm->SetBlock(i, vtu);
    }
  }

  std::vector<ttk::DiagramType> intermediateDiagrams(nDiags);
  std::vector<ttk::DiagramType> intermediateAtoms(numInputAtoms);

  double max_dimension_total = 0.0;
  for(int i = 0; i < nDiags; ++i) {

    const auto ret
      = VTUToDiagram(intermediateDiagrams[i], inputDiagrams[i], *this);

    double maxPers = this->getMaxPers(intermediateDiagrams[i]);
    double percentage = this->Percent_;

    intermediateDiagrams[i].erase(
      std::remove_if(intermediateDiagrams[i].begin(),
                     intermediateDiagrams[i].end(),
                     [maxPers, percentage](ttk::PersistencePair &t) {
                       return (t.death.sfValue - t.birth.sfValue)
                              < (percentage / 100.) * maxPers;
                     }),
      intermediateDiagrams[i].end());

    if(ret != 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
    if(max_dimension_total < maxPers) {
      max_dimension_total = maxPers;
    }
  }
  for(int i = 0; i < numInputAtoms; ++i) {
    const auto ret = VTUToDiagram(intermediateAtoms[i], inputAtoms[i], *this);
    if(ret != 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
  }

  std::vector<ttk::DiagramType> dictDiagrams;
  const int seed = this->GetSeed_();

  std::vector<std::vector<double>> vectorWeights(nDiags);
  for(size_t i = 0; i < vectorWeights.size(); ++i) {
    std::vector<double> weights(numAtom, 1. / (numAtom * 1.));
    vectorWeights[i] = std::move(weights);
  }

  std::vector<double> lossTab;
  std::vector<std::vector<double>> allLosses(nDiags);
  this->execute(intermediateDiagrams, intermediateAtoms, dictDiagrams,
                vectorWeights, seed, numAtom, lossTab, allLosses,
                this->Percent_);
  // zero-padd column name to keep Row Data columns ordered
  outputWeights->SetNumberOfRows(nDiags);

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };
  for(int i = 0; i < numAtom; ++i) {
    std::string name{"Atom"};
    zeroPad(name, numAtom, i);
    // name
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(nDiags);
    col->SetName(name.c_str());
    for(int j = 0; j < nDiags; ++j) {
      col->SetValue(j, vectorWeights[j][i]);
    }
    col->Modified();
    outputWeights->AddColumn(col);
  }

  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputDiagrams[0]->GetFieldData());
  fd->SetNumberOfTuples(nDiags);
  for(int i = 0; i < nDiags; ++i) {
    fd->SetTuple(i, 0, inputDiagrams[i]->GetFieldData());
  }

  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    outputWeights->AddColumn(fd->GetAbstractArray(i));
  }

  vtkNew<vtkFloatArray> dummy{};

  for(int i = 0; i < numAtom; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu;
    ttk::DiagramType &diagram = dictDiagrams[i];
    DiagramToVTU(vtu, diagram, dummy, *this, 3, false);
    outputDgm->SetBlock(i, vtu);
  }
  return 1;
}
