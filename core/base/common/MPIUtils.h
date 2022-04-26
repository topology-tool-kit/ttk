/// \ingroup base
/// \class ttk::MPIUtils
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
///\brief TTK base package for MPI related utilites.

#ifdef TTK_ENABLE_MPI
#include <mpi.h>
#endif

namespace ttk {
#ifdef TTK_ENABLE_MPI
  MPI_Datatype getMPIType(const float ttkNotUsed(val)) {
    return MPI_FLOAT;
  };
  MPI_Datatype getMPIType(const int ttkNotUsed(val)) {
    return MPI_INT;
  };
  MPI_Datatype getMPIType(const unsigned int ttkNotUsed(val)) {
    return MPI_UNSIGNED;
  };
  MPI_Datatype getMPIType(const double ttkNotUsed(val)) {
    return MPI_DOUBLE;
  };
  MPI_Datatype getMPIType(const long double ttkNotUsed(val)) {
    return MPI_LONG_DOUBLE;
  };
  MPI_Datatype getMPIType(const long ttkNotUsed(val)) {
    return MPI_LONG;
  };
  MPI_Datatype getMPIType(const unsigned long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG;
  };
  MPI_Datatype getMPIType(const long long ttkNotUsed(val)) {
    return MPI_LONG_LONG;
  };
  MPI_Datatype getMPIType(const unsigned long long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG_LONG;
  };

  bool isRunningWithMPI() {
    int flag_i;
    MPI_Initialized(&flag_i);
    return flag_i;
  }
#endif

} // namespace ttk
