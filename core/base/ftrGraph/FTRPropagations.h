/// \ingroup base
/// \class ttk::ftr::Propagations
/// \author charles gueunet charles.gueunet+ttk@gmail.com
/// \date 2018-07-11
///
/// \brief manage propagations for FTR Graph
///
/// \sa ttk::Triangulation
/// \sa FTRGraph.h %for a usage example.

#ifndef PROPAGATIONS_H
#define PROPAGATIONS_H

// local includes
#include "FTRAtomicVector.h"
#include "FTRDataTypes.h"
#include "FTRPropagation.h"

// c++ includes
#include <vector>

namespace ttk {
  namespace ftr {
    struct Visit {
      Propagation *prop;
      bool done;
    };

    struct Visits {
      std::vector<Visit> up, down;
    };

    // Split in one up one down ?
    class Propagations : public Allocable {
      FTRAtomicVector<Propagation *> propagations_;
      Visits visits_;

    public:
      virtual ~Propagations() {
        for(Propagation *p : propagations_) {
          delete p;
        }
      }

      void alloc() override {
        propagations_.reserve(nbElmt_);
        visits_.down.resize(nbElmt_);
        visits_.up.resize(nbElmt_);
      }

      void init() override {
        fillVector<Visit>(visits_.down, {nullptr, false});
        fillVector<Visit>(visits_.up, {nullptr, false});
      }

      // newPropagation
      // history / toVisit related function
      // localGrowth maybe :P

      // Create a new propagation starting at leaf
      Propagation *newPropagation(const idVertex leaf,
                                  const VertCompFN &comp,
                                  const bool fromMin) {
        Propagation *localProp = new Propagation(leaf, comp, fromMin);
        const auto propId = propagations_.getNext();
        propagations_[propId] = localProp;
        return localProp;
      }

      void toVisit(const idVertex v, Propagation *const prop) {
        if(prop->goUp()) {
          visits_.up[v].prop = prop;
        } else {
          visits_.down[v].prop = prop;
        }
      }

      void visit(const idVertex v, Propagation *const prop) {
        if(prop->goUp()) {
          visits_.up[v].prop = prop;
          visits_.up[v].done = true;
        } else {
          visits_.down[v].prop = prop;
          visits_.down[v].done = true;
        }
      }

      bool willVisit(const idVertex v, const Propagation *const prop) const {
        if(prop->goUp()) {
          return visits_.up[v].prop == prop;
        } else {
          return visits_.down[v].prop == prop;
        }
      }

      bool hasVisited(const idVertex v, Propagation *const prop) const {
        // return visits_[v].prop && visits_[v].prop->getId() == prop->getId()
        // && visits_[v].done; more revisit but no UF traversal. No perf impact
        // noticed
        if(prop->goUp()) {
          return visits_.up[v].prop == prop && visits_.up[v].done;
        } else {
          return visits_.down[v].prop == prop && visits_.down[v].done;
        }
      }

      // check if the opposite propagation already processed this vertex (done)
      // TODO: improve by accepting vertices planned to be visited (not the
      // .done)
      bool hasVisitedOpposite(const idVertex v, Propagation *const prop) const {
        // reversed
        bool res;
        if(prop->goUp()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
          res = visits_.down[v].done;
        } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
          res = visits_.up[v].done;
        }
        return res;
      }

      Visit visit(const idVertex v, const Propagation *const prop) const {
        if(prop->goUp()) {
          return visits_.up[v];
        } else {
          return visits_.down[v];
        }
      }

      Visit visitOpposite(const idVertex v,
                          const Propagation *const prop) const {
        if(prop->goUp()) {
          return visits_.down[v];
        } else {
          return visits_.up[v];
        }
      }
    };
  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: PROPAGATIONS_H */
