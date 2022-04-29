/// \ingroup base
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Utilities for MPI implementation.

#pragma once

#include <BaseClass.h>

#if TTK_ENABLE_MPI
#include <mpi.h>

namespace ttk {

  inline MPI_Datatype getMPIType(const float ttkNotUsed(val)) {
    return MPI_FLOAT;
  };
  inline MPI_Datatype getMPIType(const int ttkNotUsed(val)) {
    return MPI_INT;
  };
  inline MPI_Datatype getMPIType(const unsigned int ttkNotUsed(val)) {
    return MPI_UNSIGNED;
  };
  inline MPI_Datatype getMPIType(const double ttkNotUsed(val)) {
    return MPI_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long double ttkNotUsed(val)) {
    return MPI_LONG_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long ttkNotUsed(val)) {
    return MPI_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG;
  };
  inline MPI_Datatype getMPIType(const long long ttkNotUsed(val)) {
    return MPI_LONG_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG_LONG;
  };

  inline bool isRunningWithMPI() {
    int flag_i;
    MPI_Initialized(&flag_i);
    return flag_i != 0;
  }

} // namespace ttk

#endif
