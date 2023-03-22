#include <ClusteringMetrics.h>
#include <Geometry.h> // To check wheter a double is zero.
#include <cmath> // For the log2 function
#include <map>
#include <vector>

ttk::ClusteringMetrics::ClusteringMetrics() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ClusteringMetrics");
}

inline int nChoose2(int x) {
  return x * (x - 1) / 2;
}

inline int checkContingencyMatSize(const ttk::ClusteringMetrics *object,
                                   const std::vector<std::vector<int>> &matrix,
                                   const size_t nPoint) {
  if(nPoint == 0) {
    object->printErr("Error: clustering on zero points.");
    return 0;
  }

  size_t nLin = matrix.size();
  if(nLin == 0) {
    object->printErr("The provided contingency matrix is empty.\n");
    return 0;
  }
  size_t nCol = matrix[0].size();

  for(size_t i = 0; i < nLin; i++) {
    size_t curNCol = matrix[i].size();
    if(curNCol == 0) {
      object->printErr("Line " + std::to_string(i)
                       + " of the contingency matrix is empty.\n");
      return 0;
    } else if(curNCol != nCol) {
      object->printErr(
        "Line " + std::to_string(i)
        + " of the contingency matrix has wrong number of columns : "
        + std::to_string(curNCol) + " instead of " + std::to_string(nCol)
        + ".\n");
      return 0;
    }
  }
  return 1;
}

int ttk::ClusteringMetrics::computeContingencyTables(
  const int *clust1,
  const int *clust2,
  const size_t nPoint,
  std::vector<std::vector<int>> &contingencyMatrix,
  std::vector<int> &sumLin,
  std::vector<int> &sumCol) const {
  if(nPoint == 0) {
    this->printErr("Error: clustering on zero points.");
    return 0;
  }

  std::map<int, int> values1ToId, values2ToId;
  size_t nbVal1 = 0, nbVal2 = 0;
  for(size_t i = 0; i < nPoint; i++) {
    int x1 = clust1[i], x2 = clust2[i];
    auto found1 = values1ToId.find(x1), found2 = values2ToId.find(x2);

    if(found1 == values1ToId.end()) {
      values1ToId[x1] = nbVal1;
      nbVal1++;
    }

    if(found2 == values2ToId.end()) {
      values2ToId[x2] = nbVal2;
      nbVal2++;
    }
  }

  size_t nCluster1 = nbVal1, nCluster2 = nbVal2;
  contingencyMatrix.resize(nCluster1);
  for(size_t i = 0; i < nCluster1; i++)
    contingencyMatrix[i].resize(nCluster2, 0);
  sumLin.resize(nCluster1);
  sumCol.resize(nCluster2, 0);

  for(size_t i = 0; i < nPoint; i++) {
    int x1 = values1ToId[clust1[i]], x2 = values2ToId[clust2[i]];
    contingencyMatrix[x1][x2]++;
  }

  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    int sum = 0;
    for(size_t i2 = 0; i2 < nCluster2; i2++) {
      sumCol[i2] += contingencyMatrix[i1][i2];
      sum += contingencyMatrix[i1][i2];
    }

    sumLin[i1] = sum;
  }

  return 0;
}

int ttk::ClusteringMetrics::computeARI(
  const std::vector<std::vector<int>> &contingencyMatrix,
  const std::vector<int> &sumLin,
  const std::vector<int> &sumCol,
  const size_t nPoint,
  double &ariValue) const {
  if(!checkContingencyMatSize(this, contingencyMatrix, nPoint))
    return 0;

  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size();

  double sumNChooseContingency = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(+:sumNChooseContingency)
#endif // TTK_ENABLE_OPENMP
  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++)
      sumNChooseContingency += nChoose2(contingencyMatrix[i1][i2]);
  }

  double sumNChoose2_1 = 0, sumNChoose2_2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
    if(sumLin[i] == 0) {
      this->printErr("Error: the sum of a line in the contingency matrix is "
                     "zero. This should not happen.");
    }
    sumNChoose2_1 += nChoose2(sumLin[i]);
  }
  for(size_t i = 0; i < nCluster2; i++) {
    if(sumCol[i] == 0) {
      this->printErr("Error: the sum of a column in the contingency matrix is "
                     "zero. This should not happen.");
    }
    sumNChoose2_2 += nChoose2(sumCol[i]);
  }

  double numerator = sumNChooseContingency
                     - (sumNChoose2_1 * sumNChoose2_2) / nChoose2(nPoint);
  double denominator = 0.5 * (sumNChoose2_1 + sumNChoose2_2)
                       - (sumNChoose2_1 * sumNChoose2_2) / nChoose2(nPoint);
  if(denominator < ttk::Geometry::powIntTen(-DBL_DIG))
    ariValue = 1;
  else
    ariValue = numerator / denominator;

  return 0;
}

int ttk::ClusteringMetrics::computeNMI(
  const std::vector<std::vector<int>> &contingencyMatrix,
  const std::vector<int> &sumLin,
  const std::vector<int> &sumCol,
  const size_t nPoint,
  double &nmiValue) const {
  if(!checkContingencyMatSize(this, contingencyMatrix, nPoint))
    return 0;

  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size();

  double mutualInfo = 0;
  bool invalidCell = false;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(+:mutualInfo)
#endif // TTK_ENABLE_OPENMP
  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++) {
      if(contingencyMatrix[i1][i2] == 0)
        continue;
      if(sumLin[i1] == 0 || sumCol[i2] == 0) {
        this->printErr("Error: a sum of a line or a column of the contingency "
                       "matrix is zero. This should not happen.");
        invalidCell = true;
        continue;
      }

      double logArg = (double)nPoint * contingencyMatrix[i1][i2]
                      / (sumLin[i1] * sumCol[i2]);
      double curAdd = contingencyMatrix[i1][i2] * log2(logArg) / (nPoint);
      mutualInfo += curAdd;
    }
  }
  if(invalidCell)
    return 0;

  double entropy1 = 0, entropy2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
    double eltLin = (double)sumLin[i] / nPoint;
    entropy1 -= eltLin * log2(eltLin);
  }
  for(size_t i = 0; i < nCluster2; i++) {
    double eltCol = (double)sumCol[i] / nPoint;
    entropy2 -= eltCol * log2(eltCol);
  }

  nmiValue = 2 * mutualInfo / (entropy1 + entropy2);

  return 0;
}

int ttk::ClusteringMetrics::execute(const int *clustering1,
                                    const int *clustering2,
                                    const size_t n,
                                    double &nmiValue,
                                    double &ariValue) const {
  ttk::Timer timer;

  this->printMsg(ttk::debug::Separator::L1);

  std::vector<std::vector<int>> contingencyMatrix;
  std::vector<int> sumLines, sumColumns;
  computeContingencyTables(
    clustering1, clustering2, n, contingencyMatrix, sumLines, sumColumns);

  computeARI(contingencyMatrix, sumLines, sumColumns, n, ariValue);
  computeNMI(contingencyMatrix, sumLines, sumColumns, n, nmiValue);

  this->printMsg("Size of output in ttk/base = 0\n");

  this->printMsg("Computed NMI value: " + std::to_string(nmiValue) + "\n");
  this->printMsg("Computed ARI value: " + std::to_string(ariValue) + "\n");
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  return 0;
}
