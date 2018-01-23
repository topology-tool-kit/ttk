/// \ingroup base
//
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK structures for the reeb graph

#pragma once

namespace ttk
{
namespace ftr
{
   // Compute parameters (global)
   struct Params {
      bool     segm        = true;
      bool     normalize   = true;
      bool     advStats    = true;
      int      samplingLvl = 0;
   };
}
}

