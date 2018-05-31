/// \ingroup base
/// \class ttk::ComparableVector
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-06-01
///
///\brief A vector on witch == compare the range on not the size

#pragma once

#include <algorithm>
#include <vector>

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#include <typeinfo>
#endif

namespace ttk
{
   template <typename type>
   class ComparableVector : public std::vector<type>
   {
     public:
      // inherit constructors
      using std::vector<type>::vector;
      // redefine comparison operators
      bool operator==(const ComparableVector<type>& v1)
      {
         if (this->size() != v1.size()) {
            return false;
         }
         return std::equal(this->begin(), this->end(), v1.begin());
      }
   };
}
