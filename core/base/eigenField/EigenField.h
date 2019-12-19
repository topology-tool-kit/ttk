/// \ingroup base
/// \class ttk::EigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date April 2019
///
/// \brief TTK processing package for computing eigenfunctions of a
/// triangular mesh.
///
/// \sa ttkEigenField.cpp % for a usage example.

#pragma once

// base code includes
#include <Debug.h>
#include <Geometry.h>
#include <Laplacian.h>
#include <Triangulation.h>

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

namespace ttk {

  class EigenField : virtual public Debug {
  public:
    EigenField() {
      this->SetDebugMsgPrefix("EigenField");
    }

    inline void preconditionTriangulation(Triangulation *triangulation) const {
      if(triangulation != nullptr) {
        if(!triangulation->hasPreconditionedVertexNeighbors()) {
          triangulation->preconditionVertexNeighbors();
        }
        // cotan weights method needs more pre-processing
        if(!triangulation->hasPreconditionedEdgeTriangles()) {
          triangulation->preconditionEdgeTriangles();
        }
      }
    }

    template <typename T>
    int execute(Triangulation *triangulation,
                T *const outputFieldPointer,
                const unsigned int eigenNumber = 500,
                bool computeStatistics = false,
                T *const outputStatistics = nullptr) const {

      this->PrintMsg(ttk::debug::Separator::L1);

      Timer t;

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)

#ifdef TTK_ENABLE_OPENMP
      Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

      using SpMat = Eigen::SparseMatrix<T>;
      using DMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

      // number of vertices
      const auto vertexNumber = triangulation->getNumberOfVertices();

      // graph laplacian of current mesh
      SpMat lap;
      // compute graph laplacian using cotangent weights
      Laplacian::cotanWeights<T>(lap, *triangulation);
      // lap is square
      eigen_plain_assert(lap.cols() == lap.rows());

      auto n = lap.cols();
      auto m = eigenNumber;
      // threshold: minimal number of eigenpairs to get a converging solution
      const size_t minEigenNumber = 20;

      if(eigenNumber == 0) {
        // default value
        m = n / 1000;
      } else if(eigenNumber < minEigenNumber) {
        m = minEigenNumber;
      }

      Spectra::SparseSymMatProd<T> op(lap);
      Spectra::SymEigsSolver<T, Spectra::LARGEST_ALGE, decltype(op)> solver(
        &op, m, 2 * m);

      solver.init();

      // number of eigenpairs correctly computed
      int nconv = solver.compute();

      switch(solver.info()) {
        case Spectra::COMPUTATION_INFO::NUMERICAL_ISSUE:
          this->PrintMsg("Numerical Issue!", ttk::debug::Priority::ERROR);
          break;
        case Spectra::COMPUTATION_INFO::NOT_CONVERGING:
          this->PrintMsg("No Convergence! (" + std::to_string(nconv)
                           + " out of "
                           + std::to_string(eigenNumber) " values computed)",
                         ttk::debug::Priority::ERROR);
          break;
        case Spectra::COMPUTATION_INFO::NOT_COMPUTED:
          this->PrintMsg("Invalid Input!", ttk::debug::Priority::ERROR);
          break;
        default:
          break;
      }

      DMat eigenvectors = solver.eigenvectors();

      auto outputEigenFunctions = static_cast<T *>(outputFieldPointer);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(SimplexId i = 0; i < vertexNumber; ++i) {
        for(size_t j = 0; j < eigenNumber; ++j) {
          // cannot avoid copy here...
          outputEigenFunctions[i * eigenNumber + j] = eigenvectors(i, j);
        }
      }

      if(computeStatistics && outputStatistics != nullptr) {

        // number of statistics components
        const int statsComp = 4;

        auto outputStats = static_cast<T *>(outputStatistics);

        // compute statistics on eigenfunctions
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(SimplexId i = 0; i < vertexNumber; ++i) {
          auto k = i * statsComp;
          // init current tuple computation
          outputStats[k] = outputEigenFunctions[i * eigenNumber];
          outputStats[k + 1] = outputEigenFunctions[i * eigenNumber];
          outputStats[k + 2] = outputEigenFunctions[i * eigenNumber];
          outputStats[k + 3] = outputEigenFunctions[i * eigenNumber]
                               * outputEigenFunctions[i * eigenNumber];
          // loop from 1
          for(size_t j = 1; j < eigenNumber; ++j) {
            outputStats[k] = std::min<T>(
              outputEigenFunctions[i * eigenNumber + j], outputStats[k]);
            outputStats[k + 1] = std::max<T>(
              outputEigenFunctions[i * eigenNumber + j], outputStats[k]);
            outputStats[k + 2] += outputEigenFunctions[i * eigenNumber + j];
            outputStats[k + 3] += outputEigenFunctions[i * eigenNumber + j]
                                  * outputEigenFunctions[i * eigenNumber + j];
            ;
          }
        }
      }

      this->PrintMsg(ttk::debug::Separator::L2); // horizontal '-' separator
      this->PrintMsg(
        "Complete", 1, globalTimer.getElapsedTime() // global progress, time
      );

#else

      this->PrintMsg("Spectra support disabled, computation skipped!",
                     ttk::debug::Priority::WARNING);
      this->PrintMsg("Please re-compile TTK with Eigen AND Spectra support to "
                     "enable this feature.",
                     ttk::debug::Priority::WARNING);

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

      this->PrintMsg(ttk::debug::Separator::L1);

      return 0;
    }
  };

} // namespace ttk
