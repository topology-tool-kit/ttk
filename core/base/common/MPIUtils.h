/// \ingroup base
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Utilities for MPI implementation.

#pragma once

#include <BaseClass.h>
#include <vector>

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

  /**
   * @brief Gather vectors on a specific rank
   *
   * @param[out] dst Container storing vectors copied from other ranks
   * @param[in] src Data to send to \ref destRank and to store in \ref dst
   * @param[in] destRank Destination process identifier
   * @return 0 in case of success
   */
  template <typename T>
  int gatherVectors(std::vector<std::vector<T>> &dst,
                    const std::vector<T> &src,
                    const int destRank) {

    if(!ttk::isRunningWithMPI()) {
      return -1;
    }

    // src sizes (gathered on rank 0)
    std::vector<unsigned long> vecSizes{};

    if(ttk::MPIrank_ == destRank) {
      vecSizes.resize(MPIsize_);
      dst.resize(MPIsize_);
    }

    const unsigned long localSize = src.size();
    // gather src sizes on destRank
    MPI_Gather(&localSize, 1, MPI_UNSIGNED_LONG, vecSizes.data(), 1,
               MPI_UNSIGNED_LONG, destRank, MPI_COMM_WORLD);

    if(ttk::MPIrank_ == destRank) {
      // allocate dst with vecSizes
      for(int i = 0; i < ttk::MPIsize_; ++i) {
        if(i == destRank) {
          continue;
        }
        dst[i].resize(vecSizes[i]);
      }

      for(int i = 0; i < ttk::MPIsize_; ++i) {
        if(i == destRank) {
          continue;
        }
        // receive src content from other ranks
        MPI_Recv(dst[i].data(), dst[i].size(), ttk::getMPIType(src[0]), i,
                 MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      dst[destRank] = std::move(src);

    } else {
      // send src content to destRank
      MPI_Send(src.data(), src.size(), ttk::getMPIType(src[0]), destRank, 0,
               MPI_COMM_WORLD);
    }

    return 0;
  }

} // namespace ttk

#endif
