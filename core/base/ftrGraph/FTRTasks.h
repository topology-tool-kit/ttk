/// \ingroup base
/// \class ttk::ftr::Tasks
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-26
///
/// \brief TTK %FTRGraph tasks management tools
///
/// This class is a little toolset for tasks management
///
/// \sa ttk::FTRGraph

#pragma once

#include "FTRDataTypes.h"
#include <Debug.h>

#include <tuple>
#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

namespace ttk {
  namespace ftr {

    // At least one of nbTasks or grainSize need to be given
    // to use this structure
    struct TaskChunk {
      idVertex nbElemt = 0;
      idPropagation nbTasks = 0;
      idVertex grainSize = 0;

      explicit TaskChunk(const idVertex nbel) : nbElemt{nbel} {
      }
    };

    // Stateless class, toolbox
    class Tasks {
    public:
      /// This function gives the size of a chunk and the number
      /// of chunk to divide the number of element (nbElemt)
      /// in chunk of grainSize elements OR in nbTasks chunks
      /// \pre nbElemt is the number of element to visit and
      /// ONLY ONE of nbTasks or grainSize should be set.
      static std::tuple<idVertex, idPropagation>
        getChunk(const TaskChunk &params) {
#ifndef TTK_ENABLE_KAMIKAZE
        Debug dbg{};
        dbg.setDebugMsgPrefix("Task");
        if(!params.nbElemt) {
          dbg.printErr("getChunk called with nbElemnt null");
        }
        if(!params.nbTasks and !params.grainSize) {
          dbg.printErr("getChunk called with neither nbtasks nor grainSize");
        }
        if(params.nbTasks and params.grainSize) {
          dbg.printErr("getChunk called with both nbtasks and grainSize");
        }
#endif
        const idVertex grainSize = (params.grainSize)
                                     ? params.grainSize
                                     : params.nbElemt / params.nbTasks;
        const idPropagation nbTasks
          = params.nbElemt / grainSize + (params.nbElemt % grainSize != 0);

        return std::make_tuple(grainSize, nbTasks);
      }

      static idVertex getBegin(const idPropagation chunkId,
                               const idVertex grainSize,
                               const idVertex offset = 0) {
        return offset + grainSize * chunkId;
      }

      static idVertex getEnd(const idPropagation chunkId,
                             const idVertex grainSize,
                             const idVertex maxEnd = nullVertex) {
        const idVertex computedEnd = grainSize * (chunkId + 1);
        return std::min(maxEnd, computedEnd);
      }
    };

    // Explain priority explicitely
    enum PriorityLevel { Min = 0, Low, Average, High, Higher, Max };
  } // namespace ftr
} // namespace ttk
