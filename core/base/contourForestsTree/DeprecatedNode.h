/// \ingroup base
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK container representing a node of the MergeTree
//
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef NODE_H
#define NODE_H

#include <vector>

#include <Debug.h>

#include "DeprecatedDataTypes.h"
#include "ExtendedUF.h"

namespace ttk {
  namespace cf {

    class Node {
      friend class MergeTree;

    private:
      // mesh vertex where this node is
      SimplexId vertexId_;
      // For leaves, linkedNode is the saddle ending the persistance pair
      // For saddle, linked is the leaf starting the persistance pair in which
      // they are
      SimplexId linkedNode_;
      // link with superArc above and below
      std::vector<idSuperArc> vect_downSuperArcList_, vect_upSuperArcList_;
      // Won't be displayed if hidden
      bool hidden_;
      // valence down / up
      std::tuple<idSuperArc, idSuperArc> valence_;

    public:
      // -----------------
      // CONSTRUCTOR
      // -----------------
      // {

      Node(const SimplexId &id, const SimplexId &linked)
        : vertexId_(id), linkedNode_(linked), hidden_(false), valence_(0, 0) {
      }

      // }
      // -----------------
      // ACCESSOR
      // ------------------
      // {

      // Vertex id
      // ........................{

      inline SimplexId getVertexId() const {
        return vertexId_;
      }

      inline void setVertexId(const SimplexId &vertexId) {
        vertexId_ = vertexId;
      }

      // }
      // Linked node
      // ........................{

      inline const SimplexId &getOrigin(void) const {
        return linkedNode_;
      }

      inline const SimplexId &getTerminaison(void) const {
        return linkedNode_;
      }

      inline void setOrigin(const SimplexId &linked) {
        linkedNode_ = linked;
      }

      inline void setTerminaison(const SimplexId &linked) {
        linkedNode_ = linked;
      }

      // }
      // vector arcs
      // ............................{

      inline idSuperArc getNumberOfDownSuperArcs() const {
        return vect_downSuperArcList_.size();
      }

      inline idSuperArc getNumberOfUpSuperArcs() const {
        return vect_upSuperArcList_.size();
      }

      inline idSuperArc getNumberOfSuperArcs() const {
        return vect_upSuperArcList_.size() + vect_downSuperArcList_.size();
      }

      inline idSuperArc getDownSuperArcId(const idSuperArc &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(neighborId >= vect_downSuperArcList_.size()) {
          std::cerr << "[Merge Tree:Node] get down on bad neighbor !";
          std::cerr << std::endl;
          return 0;
        }
#endif
        return vect_downSuperArcList_[neighborId];
      }

      inline idSuperArc getUpSuperArcId(const idSuperArc &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(neighborId >= vect_upSuperArcList_.size()) {
          std::cerr << "[MergeTree:Node] No SuperArc to access "
                    << static_cast<unsigned>(neighborId);
          std::cerr << std::endl;
        }
#endif
        if(vect_upSuperArcList_.size() == 0) {
          return nullSuperArc;
        }
        return vect_upSuperArcList_[neighborId];
      }

      inline void addDownSuperArcId(const idSuperArc &downSuperArcId) {
        vect_downSuperArcList_.emplace_back(downSuperArcId);
      }

      inline void addUpSuperArcId(const idSuperArc &upSuperArcId) {
        vect_upSuperArcList_.emplace_back(upSuperArcId);
      }

      inline idSuperArc clearDownSuperArcs(void) {
        idSuperArc s = vect_downSuperArcList_.size();
        vect_downSuperArcList_.clear();
        return s;
      }

      inline idSuperArc clearUpSuperArcs(void) {
        idSuperArc s = vect_upSuperArcList_.size();
        vect_upSuperArcList_.clear();
        return s;
      }

      // remove the i^th arc
      inline void removeDownSuperArcPos(const idSuperArc &i) {
        vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
        vect_downSuperArcList_.pop_back();

        decDownValence();
      }

      // Find and remove the arc
      inline void removeDownSuperArc(const idSuperArc &idSa) {
        for(idSuperArc i = 0; i < vect_downSuperArcList_.size(); ++i) {
          if(vect_downSuperArcList_[i] == idSa) {
            vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
            vect_downSuperArcList_.pop_back();

            decDownValence();
            return;
          }
        }
      }

      // Find and remove the arc (better perf for young added arc)
      inline void removeDownSuperArcFromLast(const idSuperArc &idSa) {
        for(idSuperArc i = vect_downSuperArcList_.size() - 1;; --i) {
          if(vect_downSuperArcList_[i] == idSa) {
            vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
            vect_downSuperArcList_.pop_back();

            decDownValence();
            return;
          }
          if(i == 0) {
            return;
          }
        }
      }

      // Find and remove the arc
      inline void removeUpSuperArc(const idSuperArc &idSa) {
        for(idSuperArc i = 0; i < vect_upSuperArcList_.size(); ++i) {
          if(vect_upSuperArcList_[i] == idSa) {
            vect_upSuperArcList_[i] = vect_upSuperArcList_.back();
            vect_upSuperArcList_.pop_back();

            decUpValence();
            return;
          }
        }
      }

      // Find and remove the arc (better perf for young added arc)
      inline void removeUpSuperArcFromLast(const idSuperArc &idSa) {
        for(idSuperArc i = vect_upSuperArcList_.size();; --i) {
          if(vect_upSuperArcList_[i] == idSa) {
            vect_upSuperArcList_[i] = vect_upSuperArcList_.back();
            vect_upSuperArcList_.pop_back();

            decUpValence();
            return;
          }
          if(i == 0) {
            return;
          }
        }
      }

      // }
      // hidden node
      // ...........................................{

      inline bool isHidden() const {
        return hidden_;
      }

      inline bool isVisible() const {
        return !hidden_;
      }

      inline void hide() {
        hidden_ = true;
      }

      inline void setHidden(const bool local_hidden) {
        hidden_ = local_hidden;
      }

      // }
      // Valence
      // .......................................... {

      inline idSuperArc getUpValence(void) const {
        return std::get<1>(valence_);
      }

      inline idSuperArc getDownValence(void) const {
        return std::get<0>(valence_);
      }

      inline idSuperArc getValence(void) const {
        return std::get<0>(valence_) + std::get<1>(valence_);
      }

      inline void setUpValence(const idSuperArc &v) {
        std::get<1>(valence_) = v;
      }

      inline void setDownValence(const idSuperArc &v) {
        std::get<0>(valence_) = v;
      }

      inline void incUpValence(void) {
        ++std::get<1>(valence_);
      }

      inline void incDownValence(void) {
        ++std::get<0>(valence_);
      }

      inline void decUpValence(void) {
        --std::get<1>(valence_);
      }

      inline void decDownValence(void) {
        --std::get<0>(valence_);
      }

      // }

      // }
    };

  } // namespace cf
} // namespace ttk
#endif /* end of include guard: NODE_H */
