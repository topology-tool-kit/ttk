#include "DimensionReduction.h"
#include "PersistenceDiagramUtils.h"
#include <PersistenceDiagramDictionaryDecoding.h>
#include <cmath>

void ttk::PersistenceDiagramDictionaryDecoding::execute(
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  std::vector<ttk::DiagramType> &Barycenters) const {

  Timer tm{};

  std::vector<std::vector<std::vector<MatchingType>>> AllMatchingsAtoms(
    Barycenters.size());

  for(size_t i = 0; i < Barycenters.size(); ++i) {
    auto &barycenter = Barycenters[i];
    auto &weight = vectorWeights[i];
    auto &matchings = AllMatchingsAtoms[i];
    computeWeightedBarycenter(
      dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter);
  }

  this->printMsg(
    "Computed barycenters", 1.0, tm.getElapsedTime(), this->threadNumber_);
}

void ttk::PersistenceDiagramDictionaryDecoding::computeAtomsCoordinates(
  std::vector<ttk::DiagramType> &atoms,
  const std::vector<std::vector<double>> &vectorWeights,
  std::vector<std::array<double, 3>> &coords,
  std::vector<std::array<double, 3>> &trueCoords,
  std::vector<double> &xVector,
  std::vector<double> &yVector,
  std::vector<double> &zVector,
  const double spacing,
  const size_t nAtoms) const {

  if(nAtoms == 2) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(2);
    const auto distMatrix = MatrixCalculator.execute(atoms, nInputs);
    coords[0][0] = 0.;
    trueCoords[0][0] = 0.;
    coords[0][1] = 0.;
    trueCoords[0][0] = 0.;
    coords[1][0] = spacing * distMatrix[0][1];
    trueCoords[1][0] = distMatrix[0][1];
    trueCoords[1][1] = 0.;

  } else if(nAtoms == 3) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix
      = MatrixCalculator.execute(atoms, nInputs);
    coords[0][0] = 0.;
    trueCoords[0][0] = 0.;
    coords[0][1] = 0.;
    trueCoords[0][1] = 0.;
    coords[1][0] = spacing * distMatrix[0][1];
    trueCoords[1][0] = distMatrix[0][1];
    coords[1][1] = 0.;
    trueCoords[1][1] = 0.;
    double distOpposed = distMatrix[2][1];
    double firstDist = distMatrix[0][1];
    double distAdja = distMatrix[0][2];
    double alpha = std::acos(
      (distOpposed * distOpposed - firstDist * firstDist - distAdja * distAdja)
      / (-2. * firstDist * distAdja));
    coords[2][0] = spacing * distAdja * std::cos(alpha);
    trueCoords[2][0] = distAdja * std::cos(alpha);
    coords[2][1] = spacing * distAdja * std::sin(alpha);
    trueCoords[2][1] = distAdja * std::sin(alpha);

  } else if(nAtoms == 4) {
    switch(this->ProjMet) {
      case BACKEND::DICTIONARY: {
        ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
        std::array<size_t, 2> nInputs{nAtoms, 0};
        MatrixCalculator.setDos(true, true, true);
        MatrixCalculator.setThreadNumber(3);
        std::vector<std::vector<double>> distMatrix
          = MatrixCalculator.execute(atoms, nInputs);
        coords[0][0] = 0.;
        trueCoords[0][0] = 0.;
        coords[0][1] = 0.;
        trueCoords[0][1] = 0.;
        coords[1][0] = spacing * distMatrix[0][1];
        trueCoords[1][0] = distMatrix[0][1];
        coords[1][1] = 0.;
        trueCoords[1][1] = 0.;
        double distOpposed = distMatrix[2][1];
        double firstDist = distMatrix[0][1];
        double distAdja = distMatrix[0][2];
        double alpha = std::acos((distOpposed * distOpposed
                                  - firstDist * firstDist - distAdja * distAdja)
                                 / (-2. * firstDist * distAdja));
        coords[2][0] = spacing * distAdja * std::cos(alpha);
        trueCoords[2][0] = distAdja * std::cos(alpha);
        coords[2][1] = spacing * distAdja * std::sin(alpha);
        trueCoords[2][1] = distAdja * std::sin(alpha);
        double firstHeight = distMatrix[3][0];
        double secondHeight = distMatrix[3][1];
        double thirdHeight = distMatrix[3][2];
        trueCoords[3][0] = (pow(firstHeight, 2) - pow(secondHeight, 2)
                            + pow(trueCoords[1][0], 2))
                           / (2 * trueCoords[1][0]);
        trueCoords[3][1]
          = (pow(firstHeight, 2) - pow(thirdHeight, 2)
             + pow(trueCoords[2][0], 2) + pow(trueCoords[2][1], 2)
             - 2 * trueCoords[3][0] * trueCoords[2][0])
            / (2 * trueCoords[2][1]);
        trueCoords[3][2]
          = std::sqrt(pow(firstHeight, 2) - pow(trueCoords[3][0], 2)
                      - pow(trueCoords[3][1], 2));

        break;
      }

      case BACKEND::MDS: {
        ttk::DimensionReduction DimProjector;
        DimProjector.setIsInputDistanceMatrix(true);
        DimProjector.setInputNumberOfComponents(3);
        ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
        std::array<size_t, 2> nInputs{nAtoms, 0};
        MatrixCalculator.setDos(true, true, true);
        MatrixCalculator.setThreadNumber(4);
        std::vector<std::vector<double>> distMatrix
          = MatrixCalculator.execute(atoms, nInputs);
        int nRow = distMatrix.size();
        std::vector<double> matrixForProjector;
        for(int i = 0; i < nRow; ++i) {
          for(int j = 0; j < nRow; ++j) {
            matrixForProjector.push_back(distMatrix[j][i]);
          }
        }
        std::vector<std::vector<double>> coordsAtom;
        DimProjector.execute(coordsAtom, matrixForProjector, nRow, nRow);

        for(size_t i = 0; i < 3; ++i) {
          for(size_t j = 0; j < nAtoms; ++j) {
            if(i == 0) {
              trueCoords[j][0] = coordsAtom[0][j];
            } else if(i == 1) {
              trueCoords[j][1] = coordsAtom[1][j];
            } else {
              trueCoords[j][2] = coordsAtom[2][j];
            }
          }
        }

        break;
      }
    }
  } else {
    switch(this->ProjMet) {
      case BACKEND::DICTIONARY: {

        std::vector<ttk::DiagramType> dictDiagrams;
        std::vector<ttk::DiagramType> intermediateAtoms;
        std::vector<double> lossTab;
        std::vector<std::vector<double>> allLosses(nAtoms);
        const int seed = 0;
        const int m = 3;
        std::vector<std::vector<double>> tempWeights(nAtoms);
        for(size_t i = 0; i < tempWeights.size(); ++i) {
          std::vector<double> weights(m, 1. / (m * 1.));
          tempWeights[i] = std::move(weights);
        }
        ttk::PersistenceDiagramDictionary DictionaryEncoder;
        DictionaryEncoder.setUseDimReduct(false);
        DictionaryEncoder.setUseProgApproach(true);
        DictionaryEncoder.execute(atoms, atoms, dictDiagrams, tempWeights, seed,
                                  m, lossTab, allLosses, 0.);
        std::vector<std::array<double, 3>> tempCoords(3);
        std::vector<std::array<double, 3>> tempTrueCoords(3);
        ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
        std::array<size_t, 2> nInputs{3, 0};
        MatrixCalculator.setDos(true, true, true);
        MatrixCalculator.setThreadNumber(3);
        std::vector<std::vector<double>> distMatrix
          = MatrixCalculator.execute(dictDiagrams, nInputs);
        tempCoords[0][0] = 0.;
        tempTrueCoords[0][0] = 0.;
        tempCoords[0][1] = 0.;
        tempTrueCoords[0][1] = 0.;
        tempCoords[1][0] = spacing * distMatrix[0][1];
        tempTrueCoords[1][0] = distMatrix[0][1];
        tempCoords[1][1] = 0.;
        tempTrueCoords[0][1] = 0.;
        double distOpposed = distMatrix[2][1];
        double firstDist = distMatrix[0][1];
        double distAdja = distMatrix[0][2];
        double alpha = std::acos((distOpposed * distOpposed
                                  - firstDist * firstDist - distAdja * distAdja)
                                 / (-2. * firstDist * distAdja));
        tempCoords[2][0] = spacing * distAdja * std::cos(alpha);
        tempTrueCoords[2][0] = distAdja * std::cos(alpha);
        tempCoords[2][1] = spacing * distAdja * std::sin(alpha);
        tempTrueCoords[2][1] = distAdja * std::sin(alpha);
        for(int i = 0; i < 2; ++i) {
          for(size_t j = 0; j < nAtoms; ++j) {
            double temp = 0.;
            for(int iAtom = 0; iAtom < 3; ++iAtom) {
              if(i == 0) {
                temp += tempWeights[j][iAtom] * tempTrueCoords[iAtom][0];
                trueCoords[j][0] = temp;
              } else {
                temp += tempWeights[j][iAtom] * tempTrueCoords[iAtom][1];
                trueCoords[j][1] = temp;
              }
            }
          }
        }

        break;
      }

      case BACKEND::MDS: {
        ttk::DimensionReduction DimProjector;
        DimProjector.setIsInputDistanceMatrix(true);
        ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
        std::array<size_t, 2> nInputs{nAtoms, 0};
        MatrixCalculator.setDos(true, true, true);
        MatrixCalculator.setThreadNumber(3);
        std::vector<std::vector<double>> distMatrix
          = MatrixCalculator.execute(atoms, nInputs);
        int nRow = distMatrix.size();
        std::vector<double> matrixForProjector;
        for(int i = 0; i < nRow; ++i) {
          for(int j = 0; j < nRow; ++j) {
            matrixForProjector.push_back(distMatrix[j][i]);
          }
        }
        std::vector<std::vector<double>> coordsAtom;
        DimProjector.execute(coordsAtom, matrixForProjector, nRow, nRow);

        for(size_t i = 0; i < 2; ++i) {
          for(size_t j = 0; j < nAtoms; ++j) {
            if(i == 0) {
              trueCoords[j][0] = coordsAtom[0][j];
            } else {
              trueCoords[j][1] = coordsAtom[1][j];
            }
          }
        }

        break;
      }
    }
  }
  size_t nDiags = vectorWeights.size();

  for(int i = 0; i < 3; ++i) {
    for(size_t j = 0; j < nDiags; ++j) {
      double temp = 0.;
      for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        if(i == 0) {
          temp += vectorWeights[j][iAtom] * trueCoords[iAtom][0];
        } else if(i == 1) {
          temp += vectorWeights[j][iAtom] * trueCoords[iAtom][1];
        } else {
          temp += vectorWeights[j][iAtom] * trueCoords[iAtom][2];
        }
      }
      if(i == 0) {
        xVector[j] = temp;
      } else if(i == 1) {
        yVector[j] = temp;
      } else {
        zVector[j] = temp;
      }
    }
  }
}