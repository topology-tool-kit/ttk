/// \ingroup base
/// \class ttk::MergeTreeAutoencoderDecoding
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// This module defines the %MergeTreeAutoencoderDecoding class that computes
/// a decoding of merge trees or persistence diagrams given the parameters of a
/// Wasserstein Auto-Encoder.
///
/// \b Related \b publication: \n
/// "Wasserstein Auto-Encoders of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2023
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeWAE/">Merge
///   Tree Wasserstein Auto-Encoder example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramWAE/">Persistence
///   Diagram Wasserstein Auto-Encoder example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <MergeTreeAutoencoder.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The MergeTreeAutoencoderDecoding class provides methods to compute TODO
   */
  class MergeTreeAutoencoderDecoding : virtual public Debug,
                                       public MergeTreeAutoencoder {

  public:
    MergeTreeAutoencoderDecoding();

    void execute(std::vector<ttk::ftm::MergeTree<float>> &originsTrees,
                 std::vector<ttk::ftm::MergeTree<float>> &originsPrimeTrees,
                 std::vector<unsigned int *> &allRevNodeCorr,
                 std::vector<unsigned int *> &allRevNodeCorrPrime,
                 std::vector<unsigned int> &allRevNodeCorrSize,
                 std::vector<unsigned int> &allRevNodeCorrPrimeSize);

  }; // MergeTreeAutoencoderDecoding class

} // namespace ttk
