#ifndef FTRGRAPHPRINT_TEMPLATE_H
#define FTRGRAPHPRINT_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printEdge(const orderedEdge&       oEdge,
                                                  const Propagation* const localPropagation) const
      {
         std::stringstream res;
         res << "e" << std::get<2>(oEdge) << ":";
         res << "(";
         res << std::get<0>(oEdge) << " - " << std::get<1>(oEdge);
         res << ")";
         return res.str();
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printEdge(const idEdge             edgeId,
                                                  const Propagation* const localPropagation) const
      {
         return printEdge(getOrderedEdge(edgeId, localPropagation), localPropagation);
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printTriangle(const orderedTriangle&   oTriangle,
                                                      const Propagation* const localPropagation) const
      {
         std::stringstream res;
         orderedEdge       e0, e1, e2;

         e0 = getOrderedEdge(std::get<0>(oTriangle), localPropagation);
         e1 = getOrderedEdge(std::get<1>(oTriangle), localPropagation);
         e2 = getOrderedEdge(std::get<2>(oTriangle), localPropagation);

         res << "t" << std::get<3>(oTriangle) << ":";
         res << "{";
         res << " " << printEdge(e0, localPropagation);
         res << " " << printEdge(e1, localPropagation);
         res << " " << printEdge(e2, localPropagation);
         res << " }";
         return res.str();
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printTriangle(const idCell             cellId,
                                                      const Propagation* const localPropagation) const
      {
         return printTriangle(getOrderedTriangle(cellId, localPropagation), localPropagation);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::printGraph(const int verbosity) const
      {
         graph_.print(verbosity);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::printTime(DebugTimer& timer, const std::string& msg,
                                           const int lvl) const
      {
         std::ostringstream outString(std::string(lvl, ' '));
         outString << msg << timer.getElapsedTime() << std::endl;
         dMsg(std::cout, outString.str(), lvl);
      }

   }
}

#endif /* end of include guard: FTRGRAPHPRINT_TEMPLATE_H */
