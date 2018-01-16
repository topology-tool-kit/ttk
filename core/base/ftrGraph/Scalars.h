/// \ingroup base
/// \class ttk::ftr::Scalars
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-15
///
/// \brief TTK %fTRGraph scalars management
///
/// This class deal with scalar related operations: store them, compare them, ...
///
/// \sa ttk::FTRGraph

#ifndef SCALARS_H
#define SCALARS_H

// TODO implement this using AOS as it is more efficient

namespace ttk
{
   namespace ftr
   {
      class Scalars
      {
        private:
        public:
         Scalars();
         virtual ~Scalars();
      };
   }
}

#endif /* end of include guard: SCALARS_H */
