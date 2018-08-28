/// \ingroup base
/// \class ttk::ftr::SuperArc
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph arc
///
/// This class manage  arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#ifndef FTR_SUPERARC_H
#define FTR_SUPERARC_H

// local includes
#include "AtomicUF.h"
#include "DataTypesFTR.h"
#include "Scalars.h"
#include "Segmentation.h"

// c++ includes
#ifndef TTK_ENABLE_KAMIKAZE
#include<iostream>
#endif

namespace ttk
{
   namespace ftr
   {
      class Node;

      class SuperArc
      {
        private:
         idNode upNodeId_;
         idNode downNodeId_;
         AtomicUF* ufProp_;
         bool visible_;
         bool empty_;
         idSuperArc merged_;
         Segment segmentation_;
#ifndef NDEBUG
         bool fromUp_;
#endif

        public:
         SuperArc(const idNode down = nullNode, const idNode up = nullNode)
             : upNodeId_{up},
               downNodeId_{down},
               ufProp_{nullptr},
               visible_{true},
               empty_{true},
               merged_{nullSuperArc},
               segmentation_{}
#ifndef NDEBUG
              ,fromUp_{false}
#endif
         {
         }

         idNode getUpNodeId(void) const
         {
            // Caution. can be nullNode
            return upNodeId_;
         }

         void setUpNodeId(const idNode id)
         {
            upNodeId_ = id;
         }

         idNode getDownNodeId(void) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if(downNodeId_ == nullNode)
            {
               std::cerr << "[FTR Graph]: Arc have null down node" << std::endl;
            }
#endif
            return downNodeId_;
         }

         void setDownNodeId(const idNode id)
         {
            downNodeId_ = id;
         }

         Propagation* getPropagation(void) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if(!ufProp_)
            {
               std::cerr << "[FTR Graph]: Arc have null UF propagation" << std::endl;
            }
#endif
            return ufProp_->find()->getPropagation();
         }

         void setUfProp(AtomicUF* const UFprop)
         {
            ufProp_ = UFprop;
         }

         void hide(void)
         {
            visible_ = false;
         }

         bool isVisible(void) const
         {
            return visible_;
         }

         void visit(void)
         {
            empty_ = false;
         }

         bool isEmpty() const
         {
            return empty_;
         }

         void merge(const idSuperArc arc)
         {
            if (merged_ == nullSuperArc) {
               merged_ = arc;
               hide();
            }
         }

         bool merged(void) const
         {
            return merged_ != nullSuperArc;
         }

         idSuperArc mergedIn(void) const
         {
            return merged_;
         }

         void restore(void)
         {
            visible_ = true;
            merged_ = nullSuperArc;
         }

         const decltype(segmentation_)& segmentation() const
         {
            return segmentation_;
         }

         decltype(segmentation_)& segmentation()
         {
            return segmentation_;
         }

#ifndef NDEBUG
         void setFromUp(bool up)
         {
            fromUp_ = up;
         }

         bool getFromUp(void) const
         {
            return fromUp_;
         }
#endif
      };
   }
}

#endif /* end of include guard: FTR_SUPERARC_H */
