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
#include "AtomicVector.h"
#include "DataTypesFTR.h"
#include "Propagation.h"

// c++ includes
#include <vector>

namespace ttk
{
   namespace ftr
   {
      struct Visit
      {
         Propagation* prop;
         bool         done;
      };

      class Propagations : public Allocable
      {
         AtomicVector<Propagation*>  propagations_;
         std::vector<Visit>          visits_;

        public:

         virtual ~Propagations()
         {
            for (Propagation* p : propagations_) {
               delete p;
            }
         }

         void alloc() override
         {
            propagations_.reserve(nbVerts_);
            visits_.resize(nbVerts_);
         }

         void init() override
         {
            fillVector<Propagation*>(propagations_, nullptr);
            fillVector<Visit>(visits_, {nullptr, false});
         }

         // newPropagation
         // history / toVisit related function
         // localGrowth maybe :P

         // Create a new propagation starting at leaf
         Propagation* newPropagation(const idVertex leaf, VertCompFN comp, const bool fromMin)
         {
            Propagation* localProp = new Propagation(leaf, comp, fromMin);
            const auto   propId    = propagations_.getNext();
            propagations_[propId]  = localProp;
            return localProp;
         }

         void toVisit(const idVertex v, Propagation* const prop) {
            visits_[v].prop = prop;
         }

         void visit(const idVertex v, Propagation* const prop)
         {
            visits_[v].prop = prop;
            visits_[v].done = true;
         }

         bool willVisit(const idVertex v, const Propagation* const prop) const
         {
            return visits_[v].prop == prop;
         }

         bool hasVisited(const idVertex v, const Propagation* const prop) const
         {
            return visits_[v].prop == prop && visits_[v].done;
         }

         Visit visit(const idVertex v) const
         {
             return visits_[v];
         }
      };
   }
}

#endif /* end of include guard: PROPAGATIONS_H */
