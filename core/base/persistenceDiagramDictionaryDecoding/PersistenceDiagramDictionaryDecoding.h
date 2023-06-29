/// \ingroup base
/// \class ttk::PersistenceDiagramDictionaryDecoding
/// \author Keanu Sisouk <keanu.sisouk@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date Mai 2023
///
/// \brief TTK processing package for the computation of a Dictionary
/// of Persistence Diagrams and barycentric weights to approximate
/// an ensemble of Persistence Diagrams.
///
/// \b Related \b publication \n
/// "Wasserstein Dictionaries of Persistence Diagrams" \n
/// Keanu Sisouk, Julie Delon and Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa PersistenceDiagramDictionary
#pragma once

// ttk common includes
#include <Debug.h>
#include <PersistenceDiagramClustering.h>
#include <PersistenceDiagramDictionary.h>
#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>

#include <DimensionReduction.h>

namespace ttk {
  using Matrice = std::vector<std::vector<double>>;
  using VectorMatchingTuple = std::vector<MatchingType>;

  /**
   * The PersistenceDiagramDictionaryDecoding class provides methods to compute
   * for each vertex of a triangulation the average scalar value of itself and
   * its direct neighbors.
   */
  class PersistenceDiagramDictionaryDecoding : virtual public Debug {

  public:
    enum class BACKEND { MDS = 0, DICTIONARY = 1 };

    PersistenceDiagramDictionaryDecoding() {
      this->setDebugMsgPrefix("PersistenceDiagramDictionaryDecoding");
    }

    void execute(std::vector<DiagramType> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 std::vector<DiagramType> &Barycenters) const;

  protected:
    BACKEND ProjMet{BACKEND::MDS};
    bool ProgBarycenter{false};

    void computeAtomsCoordinates(
      std::vector<ttk::DiagramType> &atoms,
      const std::vector<std::vector<double>> &vectorWeights,
      std::vector<std::array<double, 3>> &coords,
      std::vector<std::array<double, 3>> &trueCoords,
      std::vector<double> &xVector,
      std::vector<double> &yVector,
      std::vector<double> &zVector,
      const double spacing,
      const size_t nAtoms) const;

  }; // PersistenceDiagramDictionaryDecoding class

} // namespace ttk
