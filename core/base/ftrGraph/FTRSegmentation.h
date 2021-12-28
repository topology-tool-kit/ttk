/// \ingroup base
/// \class ttk::ftr::Segment
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-08-02
///
///\brief TTK classes for the segmentation of an arc in the Reeb Graph
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef FTR_SEGMENTATION_H_
#define FTR_SEGMENTATION_H_

#include "FTRDataTypes.h"
#include "FTRScalars.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif
#include <sstream>

#include <list>
#include <tuple>
#include <vector>

namespace ttk {
  namespace ftr {
    using segm_it = std::vector<idVertex>::iterator;
    using segm_rev_it = std::vector<idVertex>::reverse_iterator;
    using segm_const_it = std::vector<idVertex>::const_iterator;
    using segm_const_rev_it = std::vector<idVertex>::const_reverse_iterator;

    class Segment {
    private:
      std::vector<idVertex> vertices_;

    public:
      explicit Segment(idVertex size);
      Segment();

      segm_const_it begin(void) const;
      segm_const_it end(void) const;
      segm_it begin(void);
      segm_it end(void);
      idVertex size(void) const;
      void reserve(const idVertex size);
      void emplace_back(const idVertex v);

      idVertex operator[](const size_t &idx) const;
      idVertex &operator[](const size_t &idx);
    };

  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: FTR_SEGMENTATION_H_ */
