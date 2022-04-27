/// \ingroup base
/// \class ttk::MultiresTopology
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date 2022.
///
/// \brief TTK processing package for progressive Topological Data Analysis
///
/// This package introduces a multiresolution hierarchical representation of the
/// data.
///
/// It allows the definition of efficient algorithms for approximate or
/// progressive computation in TDA.
//
/// It is used in ttk::ProgressiveTopology for the definition of
/// efficient progressive algorithms for the computation of critical points and
/// extremum-saddle persistence diagrams.
///
/// It is also used in ttk::ApproximateTopology for the approximate computation
/// of the extremum-saddle persistence diagrams with guarantees.
///
/// \b Related \b publications \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// "Fast Approximation of Persistence Diagrams with Guarantees" \n
/// Jules Vidal, Julien Tierny\n
/// IEEE Symposium on Large Data Visualization and Analysis (LDAV), 2021
///
///
/// \sa ProgressiveTopology
/// \sa ApproximateTopology

#pragma once

// base code includes
#include <DynamicTree.h>
#include <ImplicitTriangulation.h>
#include <MultiresTriangulation.h>

#include <limits>
#include <tuple>

namespace ttk {

  using triplet = std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>;
  using polarity = unsigned char;

  class MultiresTopology : public Debug {

  public:
    struct PersistencePair {
      /** first (lower) vertex id */
      ttk::SimplexId birth{};
      /** second (higher) vertex id */
      ttk::SimplexId death{};
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      ttk::SimplexId pairType{};

      PersistencePair() = default;
      PersistencePair(const SimplexId b,
                      const SimplexId d,
                      const SimplexId pType)
        : birth{b}, death{d}, pairType{pType} {
      }
    };

    MultiresTopology() {
      this->setDebugMsgPrefix("MultiresTopology");
    }

    inline void setupTriangulation(ImplicitTriangulation *const data) {
      triangulation_ = data;
      multiresTriangulation_.setTriangulation(triangulation_);
    }

    virtual void setStartingDecimationLevel(int data) {
      startingDecimationLevel_ = std::max(data, 0);
    }
    virtual void setStoppingDecimationLevel(int data) {
      stoppingDecimationLevel_ = std::max(data, 0);
    }
    virtual void setPreallocateMemory(const bool b) {
      this->preallocateMemory_ = b;
    }
    inline int getStoppingDecimationLevel() {
      return this->stoppingDecimationLevel_;
    }

    void setStartingResolutionLevel(int rl) {
      this->setStartingDecimationLevel(multiresTriangulation_.RL_to_DL(rl));
    }

    void setStoppingResolutionLevel(int rl) {
      this->setStoppingDecimationLevel(multiresTriangulation_.RL_to_DL(rl));
    }

  protected:
    // maximum link size in 3D
    static const size_t nLink_ = 27;

    using VLBoundaryType
      = std::array<std::vector<std::pair<SimplexId, SimplexId>>, nLink_>;

    void buildVertexLinkByBoundary(const SimplexId vertexId,
                                   VLBoundaryType &vlbt) const;

    char getCriticalTypeFromLink(
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link) const;

    void getValencesFromLink(
      const SimplexId vertexId,
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link,
      std::vector<polarity> &toPropageMin,
      std::vector<polarity> &toPropageMax,
      std::vector<std::vector<SimplexId>> &saddleCCMin,
      std::vector<std::vector<SimplexId>> &saddleCCMax) const;

    void getTripletsFromSaddles(
      const SimplexId vertexId,
      std::vector<triplet> &triplets,
      const std::vector<std::vector<SimplexId>> &vertexReps) const;

    void updateLinkPolarityPonctual(
      std::vector<std::pair<polarity, polarity>> &vlp) const;

    std::string resolutionInfoString();

    ImplicitTriangulation *triangulation_{};

    MultiresTriangulation multiresTriangulation_{};

    // store the two global extrema extracted from the whole dataset vertices
    // sorting operation
    mutable SimplexId globalMax_{}, globalMin_{};

    int decimationLevel_{};
    int startingDecimationLevel_{};
    int stoppingDecimationLevel_{};
    bool preallocateMemory_{true};
    std::vector<PersistencePair> CTDiagram_{};
  };
} // namespace ttk
