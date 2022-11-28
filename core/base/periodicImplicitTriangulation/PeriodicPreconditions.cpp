#include <PeriodicPreconditions.h>

int ttk::PeriodicWithPreconditions::preconditionVerticesInternal() {
  vertexCoords_.resize(vertexNumber_);

  if(dimensionality_ == 1) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexCoords_[i][0] = i;
    }
  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      vertexToPosition2d(i, p.data());
    }
  } else if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      vertexToPosition(i, p.data());
    }
  }

  return 0;
}

int ttk::PeriodicWithPreconditions::preconditionEdgesInternal() {
  edgePositions_.resize(edgeNumber_);
  edgeCoords_.resize(edgeNumber_);

  if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < esetshift_[0]) {
        edgeToPosition(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_3D;
      } else if(i < esetshift_[1]) {
        edgeToPosition(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_3D;
      } else if(i < esetshift_[2]) {
        edgeToPosition(i, 2, p.data());
        edgePositions_[i] = EdgePosition::P_3D;
      } else if(i < esetshift_[3]) {
        edgeToPosition(i, 3, p.data());
        edgePositions_[i] = EdgePosition::D1_3D;
      } else if(i < esetshift_[4]) {
        edgeToPosition(i, 4, p.data());
        edgePositions_[i] = EdgePosition::D2_3D;
      } else if(i < esetshift_[5]) {
        edgeToPosition(i, 5, p.data());
        edgePositions_[i] = EdgePosition::D3_3D;
      } else if(i < esetshift_[6]) {
        edgeToPosition(i, 6, p.data());
        edgePositions_[i] = EdgePosition::D4_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < esetshift_[0]) {
        edgeToPosition2d(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_2D;
      } else if(i < esetshift_[1]) {
        edgeToPosition2d(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_2D;
      } else if(i < esetshift_[2]) {
        edgeToPosition2d(i, 2, p.data());
        edgePositions_[i] = EdgePosition::D1_2D;
      }
    }
  }

  else if(this->dimensionality_ == 1) {
    edgePositions_[0] = EdgePosition::FIRST_EDGE_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < edgeNumber_ - 1; ++i) {
      edgePositions_[i] = EdgePosition::CENTER_1D;
    }
    edgePositions_[edgeNumber_ - 1] = EdgePosition::LAST_EDGE_1D;
  }

  edgeVertexAccelerated_.resize(edgeNumber_);

  if(isAccelerated_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = (p[1] << div_[0]) + (p[2] << div_[1]);
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = p[1] * vshift_[0] + p[2] * vshift_[1];
    }
  }

  return 0;
}

int ttk::PeriodicWithPreconditions::preconditionTrianglesInternal() {
  if(this->dimensionality_ != 3 && this->dimensionality_ != 2) {
    return 1;
  }

  trianglePositions_.resize(triangleNumber_);
  triangleCoords_.resize(triangleNumber_);

  if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      if(i < tsetshift_[0]) {
        triangleToPosition(i, 0, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::F_3D;
      } else if(i < tsetshift_[1]) {
        triangleToPosition(i, 1, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::H_3D;
      } else if(i < tsetshift_[2]) {
        triangleToPosition(i, 2, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::C_3D;
      } else if(i < tsetshift_[3]) {
        triangleToPosition(i, 3, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D1_3D;
      } else if(i < tsetshift_[4]) {
        triangleToPosition(i, 4, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D2_3D;
      } else if(i < tsetshift_[5]) {
        triangleToPosition(i, 5, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D3_3D;
      }
    }

  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleToPosition2d(i, triangleCoords_[i].data());
      if(i % 2 == 0) {
        trianglePositions_[i] = TrianglePosition::TOP_2D;
      } else {
        trianglePositions_[i] = TrianglePosition::BOTTOM_2D;
      }
    }
  }

  return 0;
}

int ttk::PeriodicWithPreconditions::preconditionTetrahedronsInternal() {

  if(this->dimensionality_ != 3) {
    return 1;
  }
  tetrahedronCoords_.resize(this->tetrahedronNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->tetrahedronNumber_; ++i) {
    this->tetrahedronToPosition(i, tetrahedronCoords_[i].data());
  }
  return 0;
}

ttk::PeriodicImplicitTriangulation::EdgePosition
  ttk::PeriodicNoPreconditions::getEdgePosition(const SimplexId e) const {

  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    if(e < esetshift_[0]) {
      edgeToPosition(e, 0, p.data());
      return EdgePosition::L_3D;
    } else if(e < esetshift_[1]) {
      edgeToPosition(e, 1, p.data());
      return EdgePosition::H_3D;
    } else if(e < esetshift_[2]) {
      edgeToPosition(e, 2, p.data());
      return EdgePosition::P_3D;
    } else if(e < esetshift_[3]) {
      edgeToPosition(e, 3, p.data());
      return EdgePosition::D1_3D;
    } else if(e < esetshift_[4]) {
      edgeToPosition(e, 4, p.data());
      return EdgePosition::D2_3D;
    } else if(e < esetshift_[5]) {
      edgeToPosition(e, 5, p.data());
      return EdgePosition::D3_3D;
    } else if(e < esetshift_[6]) {
      edgeToPosition(e, 6, p.data());
      return EdgePosition::D4_3D;
    }
  } else if(this->dimensionality_ == 2) {
    if(e < esetshift_[0]) {
      edgeToPosition2d(e, 0, p.data());
      return EdgePosition::L_2D;
    } else if(e < esetshift_[1]) {
      edgeToPosition2d(e, 1, p.data());
      return EdgePosition::H_2D;
    } else if(e < esetshift_[2]) {
      edgeToPosition2d(e, 2, p.data());
      return EdgePosition::D1_2D;
    }
  } else if(this->dimensionality_ == 1) {
    if(e == 0) {
      return EdgePosition::FIRST_EDGE_1D;
    } else if(e == this->edgeNumber_ - 1)
      return EdgePosition::CENTER_1D;
  } else {
    return EdgePosition::LAST_EDGE_1D;
  }

  return EdgePosition::CENTER_1D;
}

std::array<ttk::SimplexId, 3>
  ttk::PeriodicNoPreconditions::getEdgeCoords(const SimplexId e) const {
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition(e, 0, p.data());
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition(e, 1, p.data());
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition(e, 2, p.data());
    } else if(e < this->esetshift_[3]) {
      this->edgeToPosition(e, 3, p.data());
    } else if(e < this->esetshift_[4]) {
      this->edgeToPosition(e, 4, p.data());
    } else if(e < this->esetshift_[5]) {
      this->edgeToPosition(e, 5, p.data());
    } else if(e < this->esetshift_[6]) {
      this->edgeToPosition(e, 6, p.data());
    }

  } else if(this->dimensionality_ == 2) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition2d(e, 0, p.data());
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition2d(e, 1, p.data());
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition2d(e, 2, p.data());
    }
  }
  return p;
}

ttk::PeriodicImplicitTriangulation::TrianglePosition
  ttk::PeriodicNoPreconditions::getTrianglePosition(const SimplexId t) const {
  if(this->dimensionality_ == 2) {
    if(t % 2 == 0) {
      return TrianglePosition::TOP_2D;
    } else {
      return TrianglePosition::BOTTOM_2D;
    }
  } else if(this->dimensionality_ == 3) {
    if(t < this->tsetshift_[0]) {
      return TrianglePosition::F_3D;
    } else if(t < this->tsetshift_[1]) {
      return TrianglePosition::H_3D;
    } else if(t < this->tsetshift_[2]) {
      return TrianglePosition::C_3D;
    } else if(t < this->tsetshift_[3]) {
      return TrianglePosition::D1_3D;
    } else if(t < this->tsetshift_[4]) {
      return TrianglePosition::D2_3D;
    } else if(t < this->tsetshift_[5]) {
      return TrianglePosition::D3_3D;
    }
  }
  return TrianglePosition::C_3D;
}

std::array<ttk::SimplexId, 3>
  ttk::PeriodicNoPreconditions::getTriangleCoords(const SimplexId t) const {
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 2) {
    this->triangleToPosition2d(t, p.data());
  } else if(this->dimensionality_ == 3) {
    if(t < this->tsetshift_[0]) {
      this->triangleToPosition(t, 0, p.data());
    } else if(t < this->tsetshift_[1]) {
      this->triangleToPosition(t, 1, p.data());
    } else if(t < this->tsetshift_[2]) {
      this->triangleToPosition(t, 2, p.data());
    } else if(t < this->tsetshift_[3]) {
      this->triangleToPosition(t, 3, p.data());
    } else if(t < this->tsetshift_[4]) {
      this->triangleToPosition(t, 4, p.data());
    } else if(t < this->tsetshift_[5]) {
      this->triangleToPosition(t, 5, p.data());
    }
  }
  return p;
}
