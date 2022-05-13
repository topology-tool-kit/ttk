#include <BaseMPIClass.h>

#if TTK_ENABLE_MPI
namespace ttk {
  COMMON_EXPORTS int MPIrank_;

  BaseMPIClass::BaseMPIClass() {
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank_);
  }

} // namespace ttk
#endif