/// \ingroup base
/// \class ttk::ArrayPreconditioning
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %ArrayPreconditioning class that generates order
/// arrays from a selection of scalar field arrays.
///

#pragma once

// ttk common includes
#include <Debug.h>

#include <vector>
namespace ttk {

  /**
   * The ArrayPreconditioning class provides methods to generate order arrays
   * from a selection of scalar field arrays.
   */
  class ArrayPreconditioning : virtual public Debug {

  public:
    ArrayPreconditioning();

    template <typename DT, typename GVGID, typename GVR, typename GVLID>
    int processScalarArray(ttk::SimplexId *orderArray,
                           const DT *scalarArray,
                           const GVGID &getVertexGlobalId,
                           const GVR &getVertexRank,
                           const GVLID &getVertexLocalId,
                           const size_t nVerts,
                           const int burstSize,
                           std::vector<int> neighbors
                           = {}) const { // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(nVerts)},
      });
      this->printMsg(ttk::debug::Separator::L1);

// -----------------------------------------------------------------------
// Computing order Array
// -----------------------------------------------------------------------
#ifdef TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI()) {
        ttk::produceOrdering<DT>(orderArray, scalarArray, getVertexGlobalId,
                                 getVertexRank, getVertexLocalId, nVerts,
                                 burstSize, neighbors);
      }
#else
      this->printMsg("MPI not enabled!");
      TTK_FORCE_USE(orderArray);
      TTK_FORCE_USE(scalarArray);
      TTK_FORCE_USE(getVertexGlobalId);
      TTK_FORCE_USE(getVertexRank);
      TTK_FORCE_USE(getVertexLocalId);
      TTK_FORCE_USE(burstSize);
      TTK_FORCE_USE(neighbors);
      return 0;
#endif

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  }; // ArrayPreconditioning class

} // namespace ttk
