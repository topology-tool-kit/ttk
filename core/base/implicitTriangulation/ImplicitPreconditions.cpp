#include <ImplicitPreconditions.h>

int ttk::ImplicitWithPreconditions::preconditionVerticesInternal() {

  vertexPositions_.resize(this->vertexNumber_);
  vertexCoords_.resize(this->vertexNumber_);

  if(this->dimensionality_ == 1) {
    vertexPositions_[0] = VertexPosition::LEFT_CORNER_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < this->vertexNumber_ - 1; ++i) {
      vertexPositions_[i] = VertexPosition::CENTER_1D;
    }
    vertexPositions_[this->vertexNumber_ - 1] = VertexPosition::RIGHT_CORNER_1D;

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      this->vertexToPosition2d(i, p.data());

      if(0 < p[0] and p[0] < this->nbvoxels_[this->Di_]) {
        if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
          vertexPositions_[i] = VertexPosition::CENTER_2D;
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_EDGE_2D; // ab
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_EDGE_2D; // cd
      } else if(p[0] == 0) {
        if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
          vertexPositions_[i] = VertexPosition::LEFT_EDGE_2D; // ac
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_LEFT_CORNER_2D; // a
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_LEFT_CORNER_2D; // c
      } else {
        if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
          vertexPositions_[i] = VertexPosition::RIGHT_EDGE_2D; // bd
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_RIGHT_CORNER_2D; // b
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_RIGHT_CORNER_2D; // d
      }
    }

  } else if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      this->vertexToPosition(i, p.data());

      if(0 < p[0] and p[0] < this->nbvoxels_[0]) {
        if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::CENTER_3D;
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::FRONT_FACE_3D; // abcd
          else
            vertexPositions_[i] = VertexPosition::BACK_FACE_3D; // efgh
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_FACE_3D; // abef
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::TOP_FRONT_EDGE_3D; // ab
          else
            vertexPositions_[i] = VertexPosition::TOP_BACK_EDGE_3D; // ef
        } else {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_FACE_3D; // cdgh
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::BOTTOM_FRONT_EDGE_3D; // cd
          else
            vertexPositions_[i] = VertexPosition::BOTTOM_BACK_EDGE_3D; // gh
        }
      } else if(p[0] == 0) {
        if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::LEFT_FACE_3D; // aceg
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::LEFT_FRONT_EDGE_3D; // ac
          else
            vertexPositions_[i] = VertexPosition::LEFT_BACK_EDGE_3D; // eg
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_LEFT_EDGE_3D; // ae
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::TOP_LEFT_FRONT_CORNER_3D; // a
          else
            vertexPositions_[i] = VertexPosition::TOP_LEFT_BACK_CORNER_3D; // e
        } else {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_LEFT_EDGE_3D; // cg
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D; // c
          else
            vertexPositions_[i]
              = VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D; // g
        }
      } else {
        if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::RIGHT_FACE_3D; // bdfh
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::RIGHT_FRONT_EDGE_3D; // bd
          else
            vertexPositions_[i] = VertexPosition::RIGHT_BACK_EDGE_3D; // fh
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_RIGHT_EDGE_3D; // bf
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::TOP_RIGHT_FRONT_CORNER_3D; // b
          else
            vertexPositions_[i] = VertexPosition::TOP_RIGHT_BACK_CORNER_3D; // f
        } else {
          if(0 < p[2] and p[2] < this->nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_RIGHT_EDGE_3D; // dh
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D; // d
          else
            vertexPositions_[i]
              = VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D; // h
        }
      }
    }
  }
#ifdef TTK_ENABLE_MPI
  return this->preconditionDistributedVertices();
#endif // TTK_ENABLE_MPI

  return 0;
}

int ttk::ImplicitWithPreconditions::preconditionEdgesInternal() {
  edgePositions_.resize(this->edgeNumber_);
  edgeCoords_.resize(this->edgeNumber_);

  if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];

      if(i < this->esetshift_[0]) {
        this->edgeToPosition(i, 0, p.data());
        if(p[1] > 0 and p[1] < this->nbvoxels_[1]) {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::L_xnn_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::L_xn0_3D;
          else
            edgePositions_[i] = EdgePosition::L_xnN_3D;
        } else if(p[1] == 0) {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::L_x0n_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::L_x00_3D;
          else
            edgePositions_[i] = EdgePosition::L_x0N_3D;
        } else {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::L_xNn_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::L_xN0_3D;
          else
            edgePositions_[i] = EdgePosition::L_xNN_3D;
        }

      } else if(i < this->esetshift_[1]) {
        this->edgeToPosition(i, 1, p.data());
        if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::H_nyn_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::H_ny0_3D;
          else
            edgePositions_[i] = EdgePosition::H_nyN_3D;
        } else if(p[0] == 0) {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::H_0yn_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::H_0y0_3D;
          else
            edgePositions_[i] = EdgePosition::H_0yN_3D;
        } else {
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            edgePositions_[i] = EdgePosition::H_Nyn_3D;
          else if(p[2] == 0)
            edgePositions_[i] = EdgePosition::H_Ny0_3D;
          else
            edgePositions_[i] = EdgePosition::H_NyN_3D;
        }

      } else if(i < this->esetshift_[2]) {
        this->edgeToPosition(i, 2, p.data());
        if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
          if(p[1] > 0 and p[1] < this->nbvoxels_[1])
            edgePositions_[i] = EdgePosition::P_nnz_3D;
          else if(p[1] == 0)
            edgePositions_[i] = EdgePosition::P_n0z_3D;
          else
            edgePositions_[i] = EdgePosition::P_nNz_3D;
        } else if(p[0] == 0) {
          if(p[1] > 0 and p[1] < this->nbvoxels_[1])
            edgePositions_[i] = EdgePosition::P_0nz_3D;
          else if(p[1] == 0)
            edgePositions_[i] = EdgePosition::P_00z_3D;
          else
            edgePositions_[i] = EdgePosition::P_0Nz_3D;
        } else {
          if(p[1] > 0 and p[1] < this->nbvoxels_[1])
            edgePositions_[i] = EdgePosition::P_Nnz_3D;
          else if(p[1] == 0)
            edgePositions_[i] = EdgePosition::P_N0z_3D;
          else
            edgePositions_[i] = EdgePosition::P_NNz_3D;
        }

      } else if(i < this->esetshift_[3]) {
        this->edgeToPosition(i, 3, p.data());
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          edgePositions_[i] = EdgePosition::D1_xyn_3D;
        else if(p[2] == 0)
          edgePositions_[i] = EdgePosition::D1_xy0_3D;
        else
          edgePositions_[i] = EdgePosition::D1_xyN_3D;

      } else if(i < this->esetshift_[4]) {
        this->edgeToPosition(i, 4, p.data());
        if(p[0] > 0 and p[0] < this->nbvoxels_[0])
          edgePositions_[i] = EdgePosition::D2_nyz_3D;
        else if(p[0] == 0)
          edgePositions_[i] = EdgePosition::D2_0yz_3D;
        else
          edgePositions_[i] = EdgePosition::D2_Nyz_3D;

      } else if(i < this->esetshift_[5]) {
        this->edgeToPosition(i, 5, p.data());
        if(p[1] > 0 and p[1] < this->nbvoxels_[1])
          edgePositions_[i] = EdgePosition::D3_xnz_3D;
        else if(p[1] == 0)
          edgePositions_[i] = EdgePosition::D3_x0z_3D;
        else
          edgePositions_[i] = EdgePosition::D3_xNz_3D;

      } else if(i < this->esetshift_[6]) {
        this->edgeToPosition(i, 6, p.data());
        edgePositions_[i] = EdgePosition::D4_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];

      if(i < this->esetshift_[0]) {
        this->edgeToPosition2d(i, 0, p.data());
        if(p[1] > 0 and p[1] < this->nbvoxels_[this->Dj_])
          edgePositions_[i] = EdgePosition::L_xn_2D;
        else if(p[1] == 0)
          edgePositions_[i] = EdgePosition::L_x0_2D;
        else
          edgePositions_[i] = EdgePosition::L_xN_2D;

      } else if(i < this->esetshift_[1]) {
        this->edgeToPosition2d(i, 1, p.data());
        if(p[0] > 0 and p[0] < this->nbvoxels_[this->Di_])
          edgePositions_[i] = EdgePosition::H_ny_2D;
        else if(p[0] == 0)
          edgePositions_[i] = EdgePosition::H_0y_2D;
        else
          edgePositions_[i] = EdgePosition::H_Ny_2D;

      } else if(i < this->esetshift_[2]) {
        this->edgeToPosition2d(i, 2, p.data());
        edgePositions_[i] = EdgePosition::D1_2D;
      }
    }

  } else if(this->dimensionality_ == 1) {
    edgePositions_[0] = EdgePosition::FIRST_EDGE_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < this->edgeNumber_ - 1; ++i) {
      edgePositions_[i] = EdgePosition::CENTER_1D;
    }
    edgePositions_[this->edgeNumber_ - 1] = EdgePosition::LAST_EDGE_1D;
  }

#ifdef TTK_ENABLE_MPI
  return this->preconditionDistributedEdges();
#endif // TTK_ENABLE_MPI

  return 0;
}

int ttk::ImplicitWithPreconditions::preconditionTrianglesInternal() {

  trianglePositions_.resize(this->triangleNumber_);
  triangleCoords_.resize(this->triangleNumber_);
  if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->triangleNumber_; ++i) {
      if(i < this->tsetshift_[0]) {
        this->triangleToPosition(i, 0, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::F_3D;
      } else if(i < this->tsetshift_[1]) {
        this->triangleToPosition(i, 1, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::H_3D;
      } else if(i < this->tsetshift_[2]) {
        this->triangleToPosition(i, 2, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::C_3D;
      } else if(i < this->tsetshift_[3]) {
        this->triangleToPosition(i, 3, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D1_3D;
      } else if(i < this->tsetshift_[4]) {
        this->triangleToPosition(i, 4, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D2_3D;
      } else if(i < this->tsetshift_[5]) {
        this->triangleToPosition(i, 5, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D3_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->triangleNumber_; ++i) {
      this->triangleToPosition2d(i, triangleCoords_[i].data());
      if(i % 2 == 0) {
        trianglePositions_[i] = TrianglePosition::TOP_2D;
      } else {
        trianglePositions_[i] = TrianglePosition::BOTTOM_2D;
      }
    }
  }

#ifdef TTK_ENABLE_MPI
  return this->preconditionDistributedTriangles();
#endif // TTK_ENABLE_MPI

  return 0;
}

int ttk::ImplicitWithPreconditions::preconditionTetrahedronsInternal() {

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

ttk::ImplicitTriangulation::VertexPosition
  ttk::ImplicitNoPreconditions::getVertexPosition(const SimplexId v) const {

  if(this->dimensionality_ == 1) {
    if(v == 0) {
      return VertexPosition::LEFT_CORNER_1D;
    } else if(v == this->vertexNumber_ - 1) {
      return VertexPosition::RIGHT_CORNER_1D;
    }
    return VertexPosition::CENTER_1D;
  } else if(this->dimensionality_ == 2) {
    const auto p{this->getVertexCoords(v)};
    if(0 < p[0] and p[0] < this->nbvoxels_[this->Di_]) {
      if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
        return VertexPosition::CENTER_2D;
      else if(p[1] == 0)
        return VertexPosition::TOP_EDGE_2D; // ab
      else
        return VertexPosition::BOTTOM_EDGE_2D; // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
        return VertexPosition::LEFT_EDGE_2D; // ac
      else if(p[1] == 0)
        return VertexPosition::TOP_LEFT_CORNER_2D; // a
      else
        return VertexPosition::BOTTOM_LEFT_CORNER_2D; // c
    } else {
      if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
        return VertexPosition::RIGHT_EDGE_2D; // bd
      else if(p[1] == 0)
        return VertexPosition::TOP_RIGHT_CORNER_2D; // b
      else
        return VertexPosition::BOTTOM_RIGHT_CORNER_2D; // d
    }

  } else if(this->dimensionality_ == 3) {
    const auto p{this->getVertexCoords(v)};
    if(0 < p[0] and p[0] < this->nbvoxels_[0]) {
      if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::CENTER_3D;
        else if(p[2] == 0)
          return VertexPosition::FRONT_FACE_3D; // abcd
        else
          return VertexPosition::BACK_FACE_3D; // efgh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::TOP_FACE_3D; // abef
        else if(p[2] == 0)
          return VertexPosition::TOP_FRONT_EDGE_3D; // ab
        else
          return VertexPosition::TOP_BACK_EDGE_3D; // ef
      } else {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::BOTTOM_FACE_3D; // cdgh
        else if(p[2] == 0)
          return VertexPosition::BOTTOM_FRONT_EDGE_3D; // cd
        else
          return VertexPosition::BOTTOM_BACK_EDGE_3D; // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::LEFT_FACE_3D; // aceg
        else if(p[2] == 0)
          return VertexPosition::LEFT_FRONT_EDGE_3D; // ac
        else
          return VertexPosition::LEFT_BACK_EDGE_3D; // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::TOP_LEFT_EDGE_3D; // ae
        else if(p[2] == 0)
          return VertexPosition::TOP_LEFT_FRONT_CORNER_3D; // a
        else
          return VertexPosition::TOP_LEFT_BACK_CORNER_3D; // e
      } else {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::BOTTOM_LEFT_EDGE_3D; // cg
        else if(p[2] == 0)
          return VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D; // c
        else
          return VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D; // g
      }
    } else {
      if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::RIGHT_FACE_3D; // bdfh
        else if(p[2] == 0)
          return VertexPosition::RIGHT_FRONT_EDGE_3D; // bd
        else
          return VertexPosition::RIGHT_BACK_EDGE_3D; // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::TOP_RIGHT_EDGE_3D; // bf
        else if(p[2] == 0)
          return VertexPosition::TOP_RIGHT_FRONT_CORNER_3D; // b
        else
          return VertexPosition::TOP_RIGHT_BACK_CORNER_3D; // f
      } else {
        if(0 < p[2] and p[2] < this->nbvoxels_[2])
          return VertexPosition::BOTTOM_RIGHT_EDGE_3D; // dh
        else if(p[2] == 0)
          return VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D; // d
        else
          return VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D; // h
      }
    }
  }
  return VertexPosition::CENTER_3D;
}

ttk::ImplicitTriangulation::EdgePosition
  ttk::ImplicitNoPreconditions::getEdgePosition(const SimplexId e) const {
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {

    if(e < this->esetshift_[0]) {
      this->edgeToPosition(e, 0, p.data());
      if(p[1] > 0 and p[1] < this->nbvoxels_[1]) {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::L_xnn_3D;
        else if(p[2] == 0)
          return EdgePosition::L_xn0_3D;
        else
          return EdgePosition::L_xnN_3D;
      } else if(p[1] == 0) {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::L_x0n_3D;
        else if(p[2] == 0)
          return EdgePosition::L_x00_3D;
        else
          return EdgePosition::L_x0N_3D;
      } else {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::L_xNn_3D;
        else if(p[2] == 0)
          return EdgePosition::L_xN0_3D;
        else
          return EdgePosition::L_xNN_3D;
      }
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition(e, 1, p.data());
      if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::H_nyn_3D;
        else if(p[2] == 0)
          return EdgePosition::H_ny0_3D;
        else
          return EdgePosition::H_nyN_3D;
      } else if(p[0] == 0) {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::H_0yn_3D;
        else if(p[2] == 0)
          return EdgePosition::H_0y0_3D;
        else
          return EdgePosition::H_0yN_3D;
      } else {
        if(p[2] > 0 and p[2] < this->nbvoxels_[2])
          return EdgePosition::H_Nyn_3D;
        else if(p[2] == 0)
          return EdgePosition::H_Ny0_3D;
        else
          return EdgePosition::H_NyN_3D;
      }
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition(e, 2, p.data());
      if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
        if(p[1] > 0 and p[1] < this->nbvoxels_[1])
          return EdgePosition::P_nnz_3D;
        else if(p[1] == 0)
          return EdgePosition::P_n0z_3D;
        else
          return EdgePosition::P_nNz_3D;
      } else if(p[0] == 0) {
        if(p[1] > 0 and p[1] < this->nbvoxels_[1])
          return EdgePosition::P_0nz_3D;
        else if(p[1] == 0)
          return EdgePosition::P_00z_3D;
        else
          return EdgePosition::P_0Nz_3D;
      } else {
        if(p[1] > 0 and p[1] < this->nbvoxels_[1])
          return EdgePosition::P_Nnz_3D;
        else if(p[1] == 0)
          return EdgePosition::P_N0z_3D;
        else
          return EdgePosition::P_NNz_3D;
      }
    } else if(e < this->esetshift_[3]) {
      this->edgeToPosition(e, 3, p.data());
      if(p[2] > 0 and p[2] < this->nbvoxels_[2])
        return EdgePosition::D1_xyn_3D;
      else if(p[2] == 0)
        return EdgePosition::D1_xy0_3D;
      else
        return EdgePosition::D1_xyN_3D;
    } else if(e < this->esetshift_[4]) {
      this->edgeToPosition(e, 4, p.data());
      if(p[0] > 0 and p[0] < this->nbvoxels_[0])
        return EdgePosition::D2_nyz_3D;
      else if(p[0] == 0)
        return EdgePosition::D2_0yz_3D;
      else
        return EdgePosition::D2_Nyz_3D;
    } else if(e < this->esetshift_[5]) {
      this->edgeToPosition(e, 5, p.data());
      if(p[1] > 0 and p[1] < this->nbvoxels_[1])
        return EdgePosition::D3_xnz_3D;
      else if(p[1] == 0)
        return EdgePosition::D3_x0z_3D;
      else
        return EdgePosition::D3_xNz_3D;
    } else if(e < this->esetshift_[6]) {
      return EdgePosition::D4_3D;
    }

  } else if(this->dimensionality_ == 2) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition2d(e, 0, p.data());
      if(p[1] > 0 and p[1] < this->nbvoxels_[this->Dj_])
        return EdgePosition::L_xn_2D;
      else if(p[1] == 0)
        return EdgePosition::L_x0_2D;
      else
        return EdgePosition::L_xN_2D;
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition2d(e, 1, p.data());
      if(p[0] > 0 and p[0] < this->nbvoxels_[this->Di_])
        return EdgePosition::H_ny_2D;
      else if(p[0] == 0)
        return EdgePosition::H_0y_2D;
      else
        return EdgePosition::H_Ny_2D;
    } else if(e < this->esetshift_[2]) {
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
  ttk::ImplicitNoPreconditions::getEdgeCoords(const SimplexId e) const {
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

ttk::ImplicitTriangulation::TrianglePosition
  ttk::ImplicitNoPreconditions::getTrianglePosition(const SimplexId t) const {
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
  ttk::ImplicitNoPreconditions::getTriangleCoords(const SimplexId t) const {
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
