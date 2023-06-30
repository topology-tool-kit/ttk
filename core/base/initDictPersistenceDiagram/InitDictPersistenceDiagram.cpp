#include <InitDictPersistenceDiagram.h>

using namespace ttk;

void InitDictPersistenceDiagram::execute(
  std::vector<ttk::DiagramType> &DictDiagrams,
  const std::vector<ttk::DiagramType> &datas,
  const int nbAtoms,
  bool do_min_,
  bool do_sad_,
  bool do_max_) {

  const size_t nDiags = datas.size();

  PersistenceDiagramDistanceMatrix MatrixCalculator;
  std::array<size_t, 2> nInputs{nDiags, 0};
  MatrixCalculator.setDos(do_min_, do_sad_, do_max_);
  MatrixCalculator.setThreadNumber(this->threadNumber_);
  Matrix distMatrix = MatrixCalculator.execute(datas, nInputs);
  std::vector<double> allDistsSummed(nDiags);
  for(size_t i = 0; i < nDiags; ++i) {
    const auto &line = distMatrix[i];
    for(size_t j = 0; j < nDiags; ++j) {
      allDistsSummed[j] += line[j];
    }
  }
  std::vector<int> indices;
  int Id1 = std::max_element(allDistsSummed.begin(), allDistsSummed.end())
            - allDistsSummed.begin();
  indices.push_back(Id1);

  for(int i = 1; i < nbAtoms; ++i) {
    indices.push_back(getNextIndex(distMatrix, indices));
  }

  std::vector<double> tempDistsSummed(nbAtoms);
  for(int i = 0; i < nbAtoms; ++i) {
    tempDistsSummed[i] = allDistsSummed[indices[i]];
  }

  int FirstId = std::min_element(tempDistsSummed.begin(), tempDistsSummed.end())
                - tempDistsSummed.begin();
  const auto &temp = datas[indices[FirstId]];
  DictDiagrams.push_back(temp);
  for(int i = 0; i < nbAtoms; ++i) {
    if(i == FirstId) {
      continue;
    } else {
      const auto &temp2 = datas[indices[i]];
      DictDiagrams.push_back(temp2);
    }
  }
}

int InitDictPersistenceDiagram::getNextIndex(
  const Matrix &distMatrix, const std::vector<int> &indices) const {
  std::vector<double> allSumCumul(distMatrix.size(), 0.);
  for(size_t k = 0; k < indices.size(); ++k) {
    const auto &line = distMatrix[indices[k]];
    for(size_t i = 0; i < distMatrix.size(); ++i) {
      if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
        allSumCumul[i] += 0.;
      } else {
        allSumCumul[i] += line[i];
      }
    }
  }
  int newId = std::max_element(allSumCumul.begin(), allSumCumul.end())
              - allSumCumul.begin();
  return newId;
}
