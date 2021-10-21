/// \ingroup base
//
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK classe representing a SuperArc of a tree,
/// containing regular vertices.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef SUPERARC_H
#define SUPERARC_H

#include <list>
#include <vector>

#include <Debug.h>

#include "DeprecatedDataTypes.h"
#include "ExtendedUF.h"

namespace ttk {
  namespace cf {

    class SuperArc {
    private:
      // Extrema
      idNode downNodeId_, upNodeId_;
      // SuperArc are in charge of interCT communication
      idPartition downCT_, upCT_, replacantCT_;
      // Before stitching, an arc can cross the interface below the partition.
      // It can also cross the interface above.
      // We keep these information
      bool overlapBelow_, overlapAbove_;
      // Keep th last vertex seen by this arc
      // After the build a a merge tree, a close step is
      // done, using this field to close each root arc
      SimplexId lastVisited_;
      // Stat of this arc, if replaced...
      ComponentState state_;
      // ... use this field to know by wich other arc.
      // Caution, we do not want chained replacant !
      idSuperArc replacantId_;

      // Regular nodes in this arc
      // Vector for initialisation only (mergetree::build & simplify)
      // Use vertList and sizeVertList_ to acces Arc's vertices after combine
      std::vector<std::pair<SimplexId, bool>> vertices_;
      // initialized with vertices_.data() and vertices_.size()
      // these two variable allow to split the segmentation when
      // inserting node in the arc without moving memory
      // |------------| Arc
      // [............] Size N
      //
      // |---*--------| Arc with new node
      // [..][........] Size N1 + N2 = N
      // So we split static array with a2 = a1 + N1
      //
      // It works because vertices are sorted by construction in these arrays

      // Array version, also Retain masqued regular nodes :
      // when a node is removed from an arc of a tree,
      // we mark as masqed vertices that are in the arc we added in CT
      // and also in the staying arc to avoid duplicate.
      // _/!\_ used at combine time => create this array using sizeVertList_
      std::pair<SimplexId, bool> *vertList_;
      // The size of *vertList_
      SimplexId sizeVertList_;
#ifndef TTK_ENABLE_KAMIKAZE
      // add a size verification for global simplify step
      SimplexId allocSgm_ = -1;
#endif

    public:
      // CONSTRUCT
      // -----------------
      // {
      SuperArc(const idNode &d,
               const idNode &u,
               const bool overB,
               const bool overA,
               const unsigned char &ctd = 0,
               const unsigned char &ctu = 0,
               const size_t &resv = 0ul,
               const ComponentState &state = ComponentState::Visible)
        : downNodeId_(d), upNodeId_(u), downCT_(ctd), upCT_(ctu),
          overlapBelow_(overB), overlapAbove_(overA), lastVisited_(nullVertex),
          state_(state), replacantId_(nullSuperArc), vertList_(NULL),
          sizeVertList_(-1) {
        vertices_.reserve(resv);
      }

      // }
      // ------------------
      // ACCESSOR
      // --------------------
      // {
      //
      // node
      // .................................{

      inline const idNode &getUpNodeId(void) const {
        return upNodeId_;
      }

      inline const idNode &getDownNodeId(void) const {
        return downNodeId_;
      }

      inline void setUpNodeId(const idNode &upId) {
        upNodeId_ = upId;
      }

      inline void setDownNodeId(const idNode &downId) {
        downNodeId_ = downId;
      }

      // }
      // tree
      // .................................{
      //
      // idPartition u_char so lighter than a ref
      inline idPartition getDownCT(void) const {
        return downCT_;
      }

      inline idPartition getUpCT(void) const {
        return upCT_;
      }

      inline void setUpCT(decltype(upCT_) &ct) {
        upCT_ = ct;
      }

      // }
      // overlap
      // .................................{

      inline bool getOverlapAbove(void) const {
        return overlapAbove_;
      }

      inline bool getOverlapBelow(void) const {
        return overlapBelow_;
      }

      inline bool isCrossing(void) const {
        return overlapBelow_ || overlapAbove_;
      }

      inline void setOverlapAbove(const bool local_overlapAbove) {
        overlapAbove_ = local_overlapAbove;
      }

      inline void setOverlapBelow(const bool local_overlapBelow) {
        overlapBelow_ = local_overlapBelow;
      }

      // }
      // last vertex seen
      // .................................{

      inline const SimplexId &getLastVisited(void) const {
        return lastVisited_;
      }

      inline void setLastVisited(const SimplexId &vertId) {
        lastVisited_ = vertId;
        vertices_.emplace_back(vertId, false);
      }

      // }
      // state
      // .................................{

      inline bool isHidden(void) const {
        return state_ == ComponentState::Hidden;
      }

      inline bool isPruned(void) const {
        return state_ == ComponentState::Merged;
      }

      inline bool isMerged(void) const {
        return state_ == ComponentState::Merged;
      }

      inline bool isExternal(void) const {
        return downCT_ != upCT_;
      }

      inline bool isVisible(void) const {
        return state_ == ComponentState::Visible;
      }

      inline void hide(void) {
        state_ = ComponentState::Hidden;
      }

      inline void merge(const idSuperArc &arc, const idPartition ct = 255) {
        replacantCT_ = (ct == 255) ? upCT_ : ct;
        replacantId_ = arc;
        state_ = ComponentState::Merged;
      }

      // }
      // replacant arc/tree (merge)
      // .................................{

      inline const idSuperArc &getReplacantArcId(void) const {
        return replacantId_;
      }

      inline idPartition getReplacantCT(void) const {
        return replacantCT_;
      }

      // }
      // regular nodes (segmentation)
      // .................................{

      inline SimplexId getNumberOfRegularNodes(void) {
        return getVertSize();
      }

      inline const SimplexId &getRegularNodeId(const SimplexId &idx) {
        return getVertList()[idx].first;
      }

      inline bool isMasqued(const SimplexId &v) const {
        return vertList_[v].second;
      }

      // The std::vector

      inline SimplexId getSegmentationSize(void) const {
        return vertices_.size();
      }

      // not const for sort in simplify
      inline std::vector<std::pair<SimplexId, bool>> &getSegmentation(void) {
        return vertices_;
      }

      // The array

      inline std::pair<SimplexId, bool> *getVertList() {
        if(sizeVertList_ == -1) {
          vertList_ = vertices_.data();
          sizeVertList_ = vertices_.size();
        }
        return vertList_;
      }

      inline const SimplexId &getVertSize() {
        if(sizeVertList_ == -1) {
          vertList_ = vertices_.data();
          sizeVertList_ = vertices_.size();
        }
        return sizeVertList_;
      }

      inline void setMasqued(const SimplexId &v) {
        vertList_[v].second = true;
      }

      inline void setVertList(std::pair<SimplexId, bool> *vl) {
        vertList_ = vl;
      }

      inline void setVertSize(const SimplexId &s) {
        sizeVertList_ = s;
      }

      // append regular nodes :

      // From array : alloc should be already done
      inline void appendVertLists(
        const std::list<std::pair<SimplexId, bool> *> &vertLists,
        std::list<SimplexId> vertSizes,
        const SimplexId &totalSize) {
        // size local
        SimplexId newSize = sizeVertList_;
        if(newSize == -1)
          newSize = 0;

        // size added
        newSize += totalSize;

        // alloc
        std::pair<SimplexId, bool> *tmpVert
          = new std::pair<SimplexId, bool>[newSize];
        SimplexId pos = 0;

        // values local
        for(SimplexId i = 0; i < sizeVertList_; ++i) {
          tmpVert[pos++] = vertList_[i];
        }

        // values added
        for(std::pair<SimplexId, bool> *vertices : vertLists) {
          const SimplexId &size = vertSizes.front();
          vertSizes.pop_front();
          for(SimplexId i = 0; i < size; ++i) {
            tmpVert[pos++] = vertices[i];
          }
        }

        // set values
        vertList_ = tmpVert;
        sizeVertList_ = newSize;
      }

      // From array : alloc should be already done (chck boundary no kamikaze
      // only)
      inline int addSegmentationGlobal(const std::pair<SimplexId, bool> *arr,
                                       const SimplexId &size) {

#ifndef TTK_ENABLE_KAMIKAZE
        // cout << "size " << sizeVertList_ << " add " << size << " on " <<
        // allocSgm_ << std::endl;
        if(sizeVertList_ + size >= allocSgm_) {
          std::cerr << "SEGMENTATION SIZE PROBLEM :" << std::endl;
          std::cerr << "alloc : " << allocSgm_ << std::endl;
          std::cerr << "size : " << sizeVertList_ << std::endl;
          std::cerr << "add : " << size << std::endl;
          // gdb
          return 1;
        }
#endif

        for(SimplexId v = 0; v < size; v++) {
          if(!arr[v].second) {
            vertList_[sizeVertList_++] = arr[v];
          }
        }

        return 0;
      }

      // from std::vector
      inline void appendSegmentation(
        const std::vector<std::pair<SimplexId, bool>> &other) {
        vertices_.insert(vertices_.end(), other.begin(), other.end());
      }

      // from std::vector with a move
      inline void
        setVertices(std::vector<std::pair<SimplexId, bool>>::iterator &a,
                    std::vector<std::pair<SimplexId, bool>>::iterator b) {
        vertices_.insert(
          vertices_.end(), make_move_iterator(a), make_move_iterator(b));
      }

      // add one vertex, alloc should be already done
      inline void addSegmentationGlobal(const SimplexId &v) {
        vertList_[sizeVertList_++] = std::make_pair(v, false);
      }

      // }

      // }
      // -----------
      // RESERVE
      // ----------
      // {
      // Trick
      // During simplification we use sizeVertList to keep the number of
      // vertices that will merge in
      // materArc
      // Doing so we ensure we have only one reserve on the std::vector
      // At the end, we set sizeVertList_ back to minus one for normal
      // utilisation.

      inline void addFuturReserve(const SimplexId &nb) {
        sizeVertList_ += nb;
        // We have an offset of -1 due to the initial value of sizeVertList_
      }

      inline void makeAllocGlobal(const SimplexId &size) {
        vertList_ = new std::pair<SimplexId, bool>[size];
        sizeVertList_ = 0;
#ifndef TTK_ENABLE_KAMIKAZE
        allocSgm_ = size;
#endif
      }

      // Real reserve
      inline void makeReserve(void) {
        if(sizeVertList_ != -1) {
          // We have an offset of -1 due to the initial value of sizeVertList_
          vertices_.reserve(sizeVertList_ + 1);
          sizeVertList_ = -1;
        }
      }

      //}
    };
  } // namespace cf
} // namespace ttk

#endif /* end of include guard: SUPERARC_H */
