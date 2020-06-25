#include <DiscreteGradient.h>

using namespace std;
using namespace ttk;
using namespace dcg;

int DiscreteGradient::getDimensionality() const {
  return dimensionality_;
}

int DiscreteGradient::getNumberOfDimensions() const {
  return dimensionality_ + 1;
}

SimplexId DiscreteGradient::getNumberOfCells(const int dimension) const {
  if(dimensionality_ == 2) {
    switch(dimension) {
      case 0:
        return inputTriangulation_->getNumberOfVertices();
        break;

      case 1:
        return inputTriangulation_->getNumberOfEdges();
        break;

      case 2:
        return inputTriangulation_->getNumberOfCells();
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(dimension) {
      case 0:
        return inputTriangulation_->getNumberOfVertices();
        break;

      case 1:
        return inputTriangulation_->getNumberOfEdges();
        break;

      case 2:
        return inputTriangulation_->getNumberOfTriangles();
        break;

      case 3:
        return inputTriangulation_->getNumberOfCells();
        break;
    }
  }

  return -1;
}

inline DiscreteGradient::lowerStarType
  DiscreteGradient::lowerStar(const SimplexId a) const {
  lowerStarType res{};

  // a belongs to its lower star
  res[0].emplace_back(CellExt{0, a});

  // store lower edges
  const auto nedges = inputTriangulation_->getVertexEdgeNumber(a);
  res[1].reserve(nedges);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    inputTriangulation_->getVertexEdge(a, i, edgeId);
    SimplexId vertexId;
    inputTriangulation_->getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == a) {
      inputTriangulation_->getEdgeVertex(edgeId, 1, vertexId);
    }
    if(vertsOrder_[vertexId] < vertsOrder_[a]) {
      res[1].emplace_back(CellExt{1, edgeId, {vertexId}, {}});
    }
  }

  if(res[1].size() < 2) {
    // at least two edges in the lower star for one triangle
    return res;
  }

  const auto processTriangle
    = [&](const SimplexId triangleId, const SimplexId v0, const SimplexId v1,
          const SimplexId v2) {
        std::array<SimplexId, 3> lowVerts{};
        if(v0 == a) {
          lowVerts[0] = v1;
          lowVerts[1] = v2;
        } else if(v1 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v2;
        } else if(v2 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
        }
        if(vertsOrder_[a] > vertsOrder_[lowVerts[0]]
           && vertsOrder_[a] > vertsOrder_[lowVerts[1]]) {
          uint8_t j{}, k{};
          // store edges indices of current triangle
          std::array<uint8_t, 3> faces{};
          for(const auto &e : res[1]) {
            if(e.lowVerts_[0] == lowVerts[0] || e.lowVerts_[0] == lowVerts[1]) {
              faces[k++] = j;
            }
            j++;
          }
          res[2].emplace_back(CellExt{2, triangleId, lowVerts, faces});
        }
      };

  if(dimensionality_ == 2) {
    // store lower triangles

    // use optimised triangulation methods:
    // getVertexStar instead of getVertexTriangle
    // getCellVertex instead of getTriangleVertex
    const auto ncells = inputTriangulation_->getVertexStarNumber(a);
    res[2].reserve(ncells);
    for(SimplexId i = 0; i < ncells; ++i) {
      SimplexId cellId;
      inputTriangulation_->getVertexStar(a, i, cellId);
      SimplexId v0{}, v1{}, v2{};
      inputTriangulation_->getCellVertex(cellId, 0, v0);
      inputTriangulation_->getCellVertex(cellId, 1, v1);
      inputTriangulation_->getCellVertex(cellId, 2, v2);
      processTriangle(cellId, v0, v1, v2);
    }
  } else if(dimensionality_ == 3) {
    // store lower triangles
    const auto ntri = inputTriangulation_->getVertexTriangleNumber(a);
    res[2].reserve(ntri);
    for(SimplexId i = 0; i < ntri; i++) {
      SimplexId triangleId;
      inputTriangulation_->getVertexTriangle(a, i, triangleId);
      SimplexId v0{}, v1{}, v2{};
      inputTriangulation_->getTriangleVertex(triangleId, 0, v0);
      inputTriangulation_->getTriangleVertex(triangleId, 1, v1);
      inputTriangulation_->getTriangleVertex(triangleId, 2, v2);
      processTriangle(triangleId, v0, v1, v2);
    }

    // at least three triangles in the lower star for one tetra
    if(res[2].size() >= 3) {
      // store lower tetra
      const auto ncells = inputTriangulation_->getVertexStarNumber(a);
      res[3].reserve(ncells);
      for(SimplexId i = 0; i < ncells; ++i) {
        SimplexId cellId;
        inputTriangulation_->getVertexStar(a, i, cellId);
        std::array<SimplexId, 3> lowVerts{};
        SimplexId v0{}, v1{}, v2{}, v3{};
        inputTriangulation_->getCellVertex(cellId, 0, v0);
        inputTriangulation_->getCellVertex(cellId, 1, v1);
        inputTriangulation_->getCellVertex(cellId, 2, v2);
        inputTriangulation_->getCellVertex(cellId, 3, v3);
        if(v0 == a) {
          lowVerts[0] = v1;
          lowVerts[1] = v2;
          lowVerts[2] = v3;
        } else if(v1 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v2;
          lowVerts[2] = v3;
        } else if(v2 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
          lowVerts[2] = v3;
        } else if(v3 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
          lowVerts[2] = v2;
        }
        if(vertsOrder_[a] > vertsOrder_[lowVerts[0]]
           && vertsOrder_[a] > vertsOrder_[lowVerts[1]]
           && vertsOrder_[a] > vertsOrder_[lowVerts[2]]) {
          uint8_t j{}, k{};
          // store triangles indices of current tetra
          std::array<uint8_t, 3> faces{};
          for(const auto &t : res[2]) {
            if((t.lowVerts_[0] == lowVerts[0] || t.lowVerts_[0] == lowVerts[1]
                || t.lowVerts_[0] == lowVerts[2])
               && (t.lowVerts_[1] == lowVerts[0]
                   || t.lowVerts_[1] == lowVerts[1]
                   || t.lowVerts_[1] == lowVerts[2])) {
              faces[k++] = j;
            }
            j++;
          }

          res[3].emplace_back(CellExt{3, cellId, lowVerts, faces});
        }
      }
    }
  }

  return res;
}

inline void DiscreteGradient::pairCells(CellExt &alpha, CellExt &beta) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  char localBId{0}, localAId{0};
  SimplexId a{}, b{};

  if(beta.dim_ == 1) {

    for(SimplexId i = 0; i < 2; ++i) {
      inputTriangulation_->getEdgeVertex(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto nedges = inputTriangulation_->getVertexEdgeNumber(alpha.id_);
    for(SimplexId i = 0; i < nedges; ++i) {
      inputTriangulation_->getVertexEdge(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
      }
    }
  } else if(beta.dim_ == 2) {
    for(SimplexId i = 0; i < 3; ++i) {
      inputTriangulation_->getTriangleEdge(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto ntri = inputTriangulation_->getEdgeTriangleNumber(alpha.id_);
    for(SimplexId i = 0; i < ntri; ++i) {
      inputTriangulation_->getEdgeTriangle(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
      }
    }
  } else {
    for(SimplexId i = 0; i < 4; ++i) {
      inputTriangulation_->getCellTriangle(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto ntetra = inputTriangulation_->getTriangleStarNumber(alpha.id_);
    for(SimplexId i = 0; i < ntetra; ++i) {
      inputTriangulation_->getTriangleStar(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
      }
    }
  }
  gradient_[alpha.dim_][alpha.dim_][alpha.id_] = localBId;
  gradient_[alpha.dim_][alpha.dim_ + 1][beta.id_] = localAId;
#else
  gradient_[alpha.dim_][alpha.dim_][alpha.id_] = beta.id_;
  gradient_[alpha.dim_][alpha.dim_ + 1][beta.id_] = alpha.id_;
#endif // TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  alpha.paired_ = true;
  beta.paired_ = true;
}

int DiscreteGradient::processLowerStars() {

  /* Compute gradient */

  auto nverts = inputTriangulation_->getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId x = 0; x < nverts; x++) {

    // Comparison function for Cells inside priority queues
    const auto orderCells = [&](const CellExt &a, const CellExt &b) -> bool {
      if(a.dim_ == b.dim_) {
        // there should be a shared facet between the two cells
        // compare the vertices not in the shared facet
        if(a.dim_ == 1) {
          return vertsOrder_[a.lowVerts_[0]] > vertsOrder_[b.lowVerts_[0]];

        } else if(a.dim_ == 2) {
          const auto &m0 = a.lowVerts_[0];
          const auto &m1 = a.lowVerts_[1];
          const auto &n0 = b.lowVerts_[0];
          const auto &n1 = b.lowVerts_[1];

          if(m0 == n0) {
            return vertsOrder_[m1] > vertsOrder_[n1];
          } else if(m0 == n1) {
            return vertsOrder_[m1] > vertsOrder_[n0];
          } else if(m1 == n0) {
            return vertsOrder_[m0] > vertsOrder_[n1];
          } else if(m1 == n1) {
            return vertsOrder_[m0] > vertsOrder_[n0];
          }

        } else if(a.dim_ == 3) {
          SimplexId m{-1}, n{-1};

          const auto &m0 = a.lowVerts_[0];
          const auto &m1 = a.lowVerts_[1];
          const auto &m2 = a.lowVerts_[2];
          const auto &n0 = b.lowVerts_[0];
          const auto &n1 = b.lowVerts_[1];
          const auto &n2 = b.lowVerts_[2];

          // extract vertex of a not in b
          if(m0 != n0 && m0 != n1 && m0 != n2) {
            m = m0;
          } else if(m1 != n0 && m1 != n1 && m1 != n2) {
            m = m1;
          } else if(m2 != n0 && m2 != n1 && m2 != n2) {
            m = m2;
          }

          // extract vertex of b not in a
          if(n0 != m0 && n0 != m1 && n0 != m2) {
            n = n0;
          } else if(n1 != m0 && n1 != m1 && n1 != m2) {
            n = n1;
          } else if(n2 != m0 && n2 != m1 && n2 != m2) {
            n = n2;
          }

          return vertsOrder_[m] > vertsOrder_[n];
        }
      } else {
        // the cell of greater dimension should contain the cell of
        // smaller dimension
        return a.dim_ > b.dim_;
      }

      return false;
    };

    // Type alias for priority queues
    using pqType
      = std::priority_queue<std::reference_wrapper<CellExt>,
                            std::vector<std::reference_wrapper<CellExt>>,
                            decltype(orderCells)>;

    // Priority queues are pushed at the beginning and popped at the
    // end. To pop the minimum, elements should be sorted in a
    // decreasing order.
    pqType pqZero(orderCells), pqOne(orderCells);

    // Insert into pqOne cofacets of cell c_alpha such as numUnpairedFaces == 1
    const auto insertCofacets = [&](const CellExt &ca, lowerStarType &ls) {
      if(ca.dim_ == 1) {
        for(auto &beta : ls[2]) {
          if(ls[1][beta.faces_[0]].id_ == ca.id_
             || ls[1][beta.faces_[1]].id_ == ca.id_) {
            // edge ca belongs to triangle beta
            if(numUnpairedFacesTriangle(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }

      } else if(ca.dim_ == 2) {
        for(auto &beta : ls[3]) {
          if(ls[2][beta.faces_[0]].id_ == ca.id_
             || ls[2][beta.faces_[1]].id_ == ca.id_
             || ls[2][beta.faces_[2]].id_ == ca.id_) {
            // triangle ca belongs to tetra beta
            if(numUnpairedFacesTetra(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }
      }
    };

    auto Lx = lowerStar(x);

    // Lx[1] empty => x is a local minimum

    if(!Lx[1].empty()) {
      // get delta: 1-cell (edge) with minimal G value (steeper gradient)
      size_t minId = 0;
      for(size_t i = 1; i < Lx[1].size(); ++i) {
        const auto &a = Lx[1][minId].lowVerts_[0];
        const auto &b = Lx[1][i].lowVerts_[0];
        if(vertsOrder_[a] > vertsOrder_[b]) {
          // edge[i] < edge[0]
          minId = i;
        }
      }

      auto &c_delta = Lx[1][minId];

      // store x (0-cell) -> delta (1-cell) V-path
      pairCells(Lx[0][0], c_delta);

      // push every 1-cell in Lx that is not delta into pqZero
      for(auto &alpha : Lx[1]) {
        if(alpha.id_ != c_delta.id_) {
          pqZero.push(alpha);
        }
      }

      // push into pqOne every coface of delta in Lx (2-cells only,
      // 3-cells have not any facet paired yet) such that
      // numUnpairedFaces == 1
      insertCofacets(c_delta, Lx);

      while(!pqOne.empty() || !pqZero.empty()) {
        while(!pqOne.empty()) {
          auto &c_alpha = pqOne.top().get();
          pqOne.pop();
          auto unpairedFaces = numUnpairedFaces(c_alpha, Lx);
          if(unpairedFaces.first == 0) {
            pqZero.push(c_alpha);
          } else {
            auto &c_pair_alpha = Lx[c_alpha.dim_ - 1][unpairedFaces.second];

            // store (pair_alpha) -> (alpha) V-path
            pairCells(c_pair_alpha, c_alpha);

            // add cofaces of c_alpha and c_pair_alpha to pqOne
            insertCofacets(c_alpha, Lx);
            insertCofacets(c_pair_alpha, Lx);
          }
        }

        // skip pair_alpha from pqZero:
        // cells in pqZero are not critical if already paired
        while(!pqZero.empty() && pqZero.top().get().paired_) {
          pqZero.pop();
        }

        if(!pqZero.empty()) {
          auto &c_gamma = pqZero.top().get();
          pqZero.pop();

          // gamma is a critical cell
          // mark gamma as paired
          c_gamma.paired_ = true;

          // add cofacets of c_gamma to pqOne
          insertCofacets(c_gamma, Lx);
        }
      }
    }
  }

  return 0;
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFaces(const CellExt &c,
                                     const lowerStarType &ls) const {
  // c.dim_ cannot be <= 1
  if(c.dim_ == 2) {
    return numUnpairedFacesTriangle(c, ls);
  } else if(c.dim_ == 3) {
    return numUnpairedFacesTetra(c, ls);
  }

  return {0, -1};
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFacesTriangle(const CellExt &c,
                                             const lowerStarType &ls) const {
  // number of unpaired faces
  std::pair<size_t, SimplexId> res{0, -1};

  // loop over edge faces of triangle
  // (2 edges per triangle in lower star)
  for(size_t i = 0; i < 2; ++i) {
    if(!ls[1][c.faces_[i]].paired_) {
      res.first++;
      res.second = c.faces_[i];
    }
  }

  return res;
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFacesTetra(const CellExt &c,
                                          const lowerStarType &ls) const {
  // number of unpaired faces
  std::pair<size_t, SimplexId> res{0, -1};

  // loop over triangle faces of tetra
  for(const auto f : c.faces_) {
    if(!ls[2][f].paired_) {
      res.first++;
      res.second = f;
    }
  }

  return res;
}

CriticalType
  DiscreteGradient::criticalTypeFromCellDimension(const int dim) const {
  if(dim == 0) {
    return CriticalType::Local_minimum;
  } else if(dim == 1) {
    return CriticalType::Saddle1;
  } else if(dim == 2 && dimensionality_ == 3) {
    return CriticalType::Saddle2;
  } else if(dim == dimensionality_) {
    return CriticalType::Local_maximum;
  } else {
    return CriticalType::Regular;
  }
}

bool DiscreteGradient::isMinimum(const Cell &cell) const {
  if(cell.dim_ == 0) {
    return (gradient_[0][0][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle1(const Cell &cell) const {
  if(cell.dim_ == 1) {
    return (gradient_[0][1][cell.id_] == -1
            and gradient_[1][1][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle2(const Cell &cell) const {
  if(dimensionality_ == 3 and cell.dim_ == 2) {
    return (gradient_[1][2][cell.id_] == -1
            and gradient_[2][2][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isMaximum(const Cell &cell) const {
  if(dimensionality_ == 2 and cell.dim_ == 2) {
    return (gradient_[1][2][cell.id_] == -1);
  }

  if(dimensionality_ == 3 and cell.dim_ == 3) {
    return (gradient_[2][3][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isCellCritical(const int cellDim,
                                      const SimplexId cellId) const {
  if(dimensionality_ == 2) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1);
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1
                and gradient_[2][2][cellId] == -1);
        break;

      case 3:
        return (gradient_[2][3][cellId] == -1);
        break;
    }
  }
  return false;
}

bool DiscreteGradient::isCellCritical(const Cell &cell) const {
  return isCellCritical(cell.dim_, cell.id_);
}

bool DiscreteGradient::isBoundary(const Cell &cell) const {
  const int cellDim = cell.dim_;
  const SimplexId cellId = cell.id_;

  if(dimensionality_ == 2) {
    switch(cellDim) {
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId edgeId;
          inputTriangulation_->getCellEdge(cellId, i, edgeId);
          if(inputTriangulation_->isEdgeOnBoundary(edgeId)) {
            return true;
          }
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cellDim) {
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        return inputTriangulation_->isTriangleOnBoundary(cellId);

      case 3:
        for(int i = 0; i < 4; ++i) {
          SimplexId triangleId;
          inputTriangulation_->getCellTriangle(cellId, i, triangleId);
          if(inputTriangulation_->isTriangleOnBoundary(triangleId)) {
            return true;
          }
        }
        break;
    }
  }
  return false;
}

SimplexId DiscreteGradient::getPairedCell(const Cell &cell,
                                          bool isReverse) const {
  SimplexId id{-1};
  if(dimensionality_ == 2) {
    switch(cell.dim_) {
      case 0:
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getVertexEdge(
          cell.id_, gradient_[0][0][cell.id_], id);
#else
        return gradient_[0][0][cell.id_];
#endif
        break;

      case 1:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getEdgeVertex(
            cell.id_, gradient_[0][1][cell.id_], id);
          return id;
#else
          return gradient_[0][1][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getEdgeStar(
          cell.id_, gradient_[1][1][cell.id_], id);
#else
        return gradient_[1][1][cell.id_];
#endif
        break;

      case 2:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getCellEdge(
            cell.id_, gradient_[1][2][cell.id_], id);
          return id;
#else
          return gradient_[1][2][cell.id_];
#endif
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cell.dim_) {
      case 0:
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getVertexEdge(
          cell.id_, gradient_[0][0][cell.id_], id);
#else
        return gradient_[0][0][cell.id_];
#endif
        break;

      case 1:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getEdgeVertex(
            cell.id_, gradient_[0][1][cell.id_], id);
          return id;
#else
          return gradient_[0][1][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getEdgeTriangle(
          cell.id_, gradient_[1][1][cell.id_], id);
#else
        return gradient_[1][1][cell.id_];
#endif
        break;

      case 2:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getTriangleEdge(
            cell.id_, gradient_[1][2][cell.id_], id);
          return id;
#else
          return gradient_[1][2][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getTriangleStar(
          cell.id_, gradient_[2][2][cell.id_], id);
#else
        return gradient_[2][2][cell.id_];
#endif
        break;

      case 3:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getCellTriangle(
            cell.id_, gradient_[2][3][cell.id_], id);
          return id;
#else
          return gradient_[2][3][cell.id_];
#endif
        }
        break;
    }
  }

  return id;
}

int DiscreteGradient::getCriticalPoints(vector<Cell> &criticalPoints) const {
  // foreach dimension
  const int numberOfDimensions = getNumberOfDimensions();
  for(int i = 0; i < numberOfDimensions; ++i) {

    // foreach cell of that dimension
    const SimplexId numberOfCells = getNumberOfCells(i);
    for(SimplexId j = 0; j < numberOfCells; ++j) {
      const Cell cell(i, j);

      if(isCellCritical(cell)) {
        criticalPoints.push_back(cell);
      }
    }
  }

  return 0;
}

int DiscreteGradient::getDescendingPath(const Cell &cell,
                                        vector<Cell> &vpath) const {
  if(dimensionality_ == 2) {
    if(cell.dim_ == 0) {
      // assume that cellId is a vertex
      SimplexId currentId = cell.id_;
      SimplexId connectedEdgeId;
      do {
        // add a vertex
        const Cell vertex(0, currentId);
        vpath.push_back(vertex);

        if(isCellCritical(vertex)) {
          break;
        }

        connectedEdgeId = getPairedCell(vertex);
        if(connectedEdgeId == -1) {
          break;
        }

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge)) {
          break;
        }

        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId != currentId) {
            currentId = vertexId;
            break;
          }
        }

      } while(connectedEdgeId != -1);
    }
  } else if(dimensionality_ == 3) {
    if(cell.dim_ == 0) {
      // assume that cellId is a vertex
      SimplexId currentId = cell.id_;
      SimplexId connectedEdgeId;
      do {
        // add a vertex
        const Cell vertex(0, currentId);
        vpath.push_back(vertex);

        if(isCellCritical(vertex)) {
          break;
        }

        connectedEdgeId = getPairedCell(vertex);
        if(connectedEdgeId == -1) {
          break;
        }

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge)) {
          break;
        }

        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId != currentId) {
            currentId = vertexId;
            break;
          }
        }

      } while(connectedEdgeId != -1);
    }
  }

  return 0;
}

bool DiscreteGradient::getDescendingPathThroughWall(
  const Cell &saddle2,
  const Cell &saddle1,
  const vector<bool> &isVisited,
  vector<Cell> *const vpath,
  const bool stopIfMultiConnected,
  const bool enableCycleDetector) const {
  // debug
  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfEdges, 0);
  }

  if(dimensionality_ == 3) {
    // add the 2-saddle to the path
    if(vpath != nullptr) {
      vpath->push_back(saddle2);
    }

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        inputTriangulation_->getTriangleEdge(saddle2.id_, i, edgeId);
        if(isVisited[edgeId]) {
          // saddle2 can be adjacent to saddle1 on the wall
          if(isSaddle1(Cell(1, edgeId))) {
            if(vpath != nullptr) {
              vpath->push_back(Cell(1, edgeId));
            }
            return 0;
          }

          currentId = edgeId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }
    }

    int oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
          cout << "[DiscreteGradient] Error : cycle detected on the wall of "
                  "1-saddleId="
               << saddle1.id_ << endl;
          break;
        }
      }

      oldId = currentId;

      // add an edge
      const Cell edge(1, currentId);
      if(vpath != nullptr) {
        vpath->push_back(edge);
      }

      if(isCellCritical(edge)) {
        break;
      }

      const SimplexId connectedTriangleId = getPairedCell(edge);

      // add a triangle
      const Cell triangle(2, connectedTriangleId);
      if(vpath != nullptr) {
        vpath->push_back(triangle);
      }

      if(isCellCritical(triangle)) {
        break;
      }

      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        inputTriangulation_->getTriangleEdge(connectedTriangleId, i, edgeId);

        if(isVisited[edgeId] and edgeId != oldId) {
          currentId = edgeId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return false;
}

int DiscreteGradient::getAscendingPath(const Cell &cell,
                                       vector<Cell> &vpath,
                                       const bool enableCycleDetector) const {

  const SimplexId numberOfCells = inputTriangulation_->getNumberOfCells();
  vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfCells, 0);
  }

  if(dimensionality_ == 2) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {
        oldId = currentId;

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) {
          break;
        }

        const SimplexId connectedEdgeId = getPairedCell(triangle, true);
        if(connectedEdgeId == -1) {
          break;
        }

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge)) {
          break;
        }

        const SimplexId starNumber
          = inputTriangulation_->getEdgeStarNumber(connectedEdgeId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          inputTriangulation_->getEdgeStar(connectedEdgeId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  } else if(dimensionality_ == 3) {
    if(cell.dim_ == 3) {
      // assume that cellId is a tetra
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {

        // debug
        if(enableCycleDetector) {
          if(isCycle[currentId] == 0) {
            isCycle[currentId] = 1;
          } else {
            cout << "[DiscreteGradient] Error : cycle detected in the path "
                    "from tetraId="
                 << cell.id_ << endl;
            break;
          }
        }

        oldId = currentId;

        // add a tetra
        const Cell tetra(3, currentId);
        vpath.push_back(tetra);

        if(isCellCritical(tetra)) {
          break;
        }

        const SimplexId connectedTriangleId = getPairedCell(tetra, true);
        if(connectedTriangleId == -1) {
          break;
        }

        // add a triangle
        const Cell triangle(2, connectedTriangleId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) {
          break;
        }

        const SimplexId starNumber
          = inputTriangulation_->getTriangleStarNumber(connectedTriangleId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          inputTriangulation_->getTriangleStar(connectedTriangleId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  }

  return 0;
}

bool DiscreteGradient::getAscendingPathThroughWall(
  const Cell &saddle1,
  const Cell &saddle2,
  const vector<bool> &isVisited,
  vector<Cell> *const vpath,
  const bool stopIfMultiConnected,
  const bool enableCycleDetector) const {
  // debug
  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
  vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfTriangles, 0);
  }

  if(dimensionality_ == 3) {
    // add the 1-saddle to the path
    if(vpath != nullptr) {
      vpath->push_back(saddle1);
    }

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      const SimplexId triangleNumber
        = inputTriangulation_->getEdgeTriangleNumber(saddle1.id_);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        inputTriangulation_->getEdgeTriangle(saddle1.id_, i, triangleId);
        if(isVisited[triangleId]) {
          // saddle1 can be adjacent to saddle2 on the wall
          if(isSaddle2(Cell(2, triangleId))) {
            if(vpath != nullptr) {
              vpath->push_back(Cell(2, triangleId));
            }
            return false;
          }

          currentId = triangleId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }
    }

    SimplexId oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
          cout << "[DiscreteGradient] Error : cycle detected on the wall of "
                  "2-saddleId="
               << saddle2.id_ << endl;
          break;
        }
      }

      oldId = currentId;

      // add a triangle
      const Cell triangle(2, currentId);
      if(vpath != nullptr) {
        vpath->push_back(triangle);
      }

      if(isCellCritical(triangle)) {
        break;
      }

      const SimplexId connectedEdgeId = getPairedCell(triangle, true);

      // add an edge
      const Cell edge(1, connectedEdgeId);
      if(vpath != nullptr) {
        vpath->push_back(edge);
      }

      if(isCellCritical(edge)) {
        break;
      }

      int nconnections = 0;
      const SimplexId triangleNumber
        = inputTriangulation_->getEdgeTriangleNumber(connectedEdgeId);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        inputTriangulation_->getEdgeTriangle(connectedEdgeId, i, triangleId);

        if(isVisited[triangleId] and triangleId != oldId) {
          currentId = triangleId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return false;
}

int DiscreteGradient::getDescendingWall(const Cell &cell,
                                        VisitedMask &mask,
                                        vector<Cell> *const wall,
                                        set<SimplexId> *const saddles) const {
  if(dimensionality_ == 3) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      const SimplexId originId = cell.id_;

      queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId triangleId = bfs.front();
        bfs.pop();

        if(!mask.isVisited_[triangleId]) {
          mask.isVisited_[triangleId] = true;
          mask.visitedIds_.emplace_back(triangleId);

          // add the triangle
          if(wall != nullptr) {
            wall->push_back(Cell(2, triangleId));
          }

          for(int j = 0; j < 3; ++j) {
            SimplexId edgeId;
            inputTriangulation_->getTriangleEdge(triangleId, j, edgeId);

            if((saddles != nullptr) and isSaddle1(Cell(1, edgeId))) {
              saddles->insert(edgeId);
            }

            const SimplexId pairedCellId = getPairedCell(Cell(1, edgeId));

            if(pairedCellId != -1 and pairedCellId != triangleId) {
              bfs.push(pairedCellId);
            }
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::getAscendingWall(const Cell &cell,
                                       VisitedMask &mask,
                                       vector<Cell> *const wall,
                                       set<SimplexId> *const saddles) const {
  if(dimensionality_ == 3) {
    if(cell.dim_ == 1) {
      // assume that cellId is an edge
      const SimplexId originId = cell.id_;

      queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId edgeId = bfs.front();
        bfs.pop();

        if(!mask.isVisited_[edgeId]) {
          mask.isVisited_[edgeId] = true;
          mask.visitedIds_.emplace_back(edgeId);

          // add the edge
          if(wall != nullptr) {
            wall->push_back(Cell(1, edgeId));
          }

          const SimplexId triangleNumber
            = inputTriangulation_->getEdgeTriangleNumber(edgeId);
          for(SimplexId j = 0; j < triangleNumber; ++j) {
            SimplexId triangleId;
            inputTriangulation_->getEdgeTriangle(edgeId, j, triangleId);

            if((saddles != nullptr) and isSaddle2(Cell(2, triangleId))) {
              saddles->insert(triangleId);
            }

            const SimplexId pairedCellId
              = getPairedCell(Cell(2, triangleId), true);

            if(pairedCellId != -1 and pairedCellId != edgeId) {
              bfs.push(pairedCellId);
            }
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPath(const vector<Cell> &vpath) {
  if(dimensionality_ == 2) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getCellEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeStarNumber(edgeId); ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeStar(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][2][triangleId] = edgeId;
      gradient_[1][1][edgeId] = triangleId;
#endif
    }
  } else if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId tetraId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 4; ++k) {
        SimplexId tmp;
        inputTriangulation_->getCellTriangle(tetraId, k, tmp);
        if(tmp == triangleId) {
          gradient_[2][3][tetraId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getTriangleStarNumber(triangleId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleStar(triangleId, k, tmp);
        if(tmp == tetraId) {
          gradient_[2][2][triangleId] = k;
          break;
        }
      }
#else
      gradient_[2][3][tetraId] = triangleId;
      gradient_[2][2][triangleId] = tetraId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPathOnWall(const vector<Cell> &vpath) {
  if(dimensionality_ == 3) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeTriangleNumber(edgeId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][2][triangleId] = edgeId;
      gradient_[1][1][edgeId] = triangleId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::reverseDescendingPathOnWall(const vector<Cell> &vpath) {
  if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId edgeId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeTriangleNumber(edgeId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][1][edgeId] = triangleId;
      gradient_[1][2][triangleId] = edgeId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::getEdgeIncenter(const SimplexId edgeId,
                                      float incenter[3]) const {
  return inputTriangulation_->getEdgeIncenter(edgeId, incenter);
}

int DiscreteGradient::getTriangleIncenter(const SimplexId triangleId,
                                          float incenter[3]) const {
  return inputTriangulation_->getTriangleIncenter(triangleId, incenter);
}

int DiscreteGradient::getTetraIncenter(const SimplexId tetraId,
                                       float incenter[3]) const {
  return inputTriangulation_->getTetraIncenter(tetraId, incenter);
}

int DiscreteGradient::getCriticalPointMap(
  const vector<pair<SimplexId, char>> &criticalPoints, vector<char> &isPL) {
  isPL.resize(numberOfVertices_);
  std::fill(isPL.begin(), isPL.end(), 0);
  for(pair<SimplexId, char> criticalPoint : criticalPoints) {
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    isPL[criticalPointId] = criticalPointType;
  }

  return 0;
}

ttk::SimplexId DiscreteGradient::getCellGreaterVertex(const Cell c) const {

  auto cellDim = c.dim_;
  auto cellId = c.id_;

  SimplexId vertexId = -1;
  if(cellDim == 0) {
    vertexId = cellId;
  }

  else if(cellDim == 1) {
    SimplexId v0;
    SimplexId v1;
    inputTriangulation_->getEdgeVertex(cellId, 0, v0);
    inputTriangulation_->getEdgeVertex(cellId, 1, v1);

    if(vertsOrder_[v0] > vertsOrder_[v1]) {
      vertexId = v0;
    } else {
      vertexId = v1;
    }
  }

  else if(cellDim == 2) {
    SimplexId v0{}, v1{}, v2{};
    inputTriangulation_->getTriangleVertex(cellId, 0, v0);
    inputTriangulation_->getTriangleVertex(cellId, 1, v1);
    inputTriangulation_->getTriangleVertex(cellId, 2, v2);
    if(vertsOrder_[v0] > vertsOrder_[v1] && vertsOrder_[v0] > vertsOrder_[v2]) {
      vertexId = v0;
    } else if(vertsOrder_[v1] > vertsOrder_[v0]
              && vertsOrder_[v1] > vertsOrder_[v2]) {
      vertexId = v1;
    } else {
      vertexId = v2;
    }
  }

  else if(cellDim == 3) {
    SimplexId v0{}, v1{}, v2{}, v3{};
    inputTriangulation_->getCellVertex(cellId, 0, v0);
    inputTriangulation_->getCellVertex(cellId, 1, v1);
    inputTriangulation_->getCellVertex(cellId, 2, v2);
    inputTriangulation_->getCellVertex(cellId, 3, v3);
    if(vertsOrder_[v0] > vertsOrder_[v1] && vertsOrder_[v0] > vertsOrder_[v2]
       && vertsOrder_[v0] > vertsOrder_[v3]) {
      vertexId = v0;
    } else if(vertsOrder_[v1] > vertsOrder_[v0]
              && vertsOrder_[v1] > vertsOrder_[v2]
              && vertsOrder_[v1] > vertsOrder_[v3]) {
      vertexId = v1;
    } else if(vertsOrder_[v2] > vertsOrder_[v0]
              && vertsOrder_[v2] > vertsOrder_[v1]
              && vertsOrder_[v2] > vertsOrder_[v3]) {
      vertexId = v2;
    } else {
      vertexId = v3;
    }
  }
  return vertexId;
}

int DiscreteGradient::setManifoldSize(
  const std::vector<Cell> &criticalPoints,
  const std::vector<size_t> &nCriticalPointsByDim,
  const std::vector<SimplexId> &maxSeeds,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold) const {

  const auto nCritPoints = criticalPoints.size();
  const auto nDimensions = getNumberOfDimensions();

  if(outputCriticalPoints_points_manifoldSize_ == nullptr) {
    return 1;
  }

  outputCriticalPoints_points_manifoldSize_->resize(nCritPoints, 0);

  // pre-compute size of descending manifold cells
  std::map<SimplexId, size_t> descendingCellsSize{};
  for(SimplexId i = 0; i < numberOfVertices_; ++i) {
    descendingCellsSize[descendingManifold[i]]++;
  }

  // minima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCriticalPointsByDim[0]; ++i) {
    const Cell &cell = criticalPoints[i];
    const SimplexId seedId = descendingManifold[cell.id_];
    const SimplexId manifoldSize = descendingCellsSize[seedId];
    (*outputCriticalPoints_points_manifoldSize_)[i] = manifoldSize;
  }

  // index of first maximum in critical points array
  size_t nFirstMaximum{};
  for(int i = 0; i < nDimensions - 1; ++i) {
    nFirstMaximum += nCriticalPointsByDim[i];
  }

  // pre-compute size of ascending manifold cells
  std::map<SimplexId, size_t> ascendingCellsSize{};
  for(SimplexId i = 0; i < numberOfVertices_; ++i) {
    ascendingCellsSize[ascendingManifold[i]]++;
  }

  // pre-compute maximum SimplexId -> index in maxSeeds
  std::map<SimplexId, size_t> seedsPos{};
  for(size_t i = 0; i < maxSeeds.size(); ++i) {
    seedsPos[maxSeeds[i]] = i;
  }

  // maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = nFirstMaximum; i < nCritPoints; ++i) {
    const Cell &cell = criticalPoints[i];
    if(seedsPos.find(cell.id_) != seedsPos.end()) {
      const auto seedId = seedsPos[cell.id_];
      const SimplexId manifoldSize = ascendingCellsSize[seedId];
      (*outputCriticalPoints_points_manifoldSize_)[i] = manifoldSize;
    }
  }

  return 0;
}

int DiscreteGradient::setGradientGlyphs(
  SimplexId &numberOfPoints,
  std::vector<float> &points,
  std::vector<char> &points_pairOrigins,
  SimplexId &numberOfCells,
  std::vector<SimplexId> &cells,
  std::vector<char> &cells_pairTypes) const {

  SimplexId pointId{};
  SimplexId cellId{};

  // foreach dimension
  const int nDimensions = getNumberOfDimensions();
  for(int i = 0; i < nDimensions - 1; ++i) {
    // foreach cell of that dimension
    const SimplexId nCells = getNumberOfCells(i);
    for(SimplexId j = 0; j < nCells; ++j) {
      const Cell cell(i, j);

      const SimplexId pairedCellId = getPairedCell(cell);
      if(pairedCellId != -1) {
        // get gradient pair
        const int pairedCellDim = i + 1;

        const Cell pairedCell(pairedCellDim, pairedCellId);

        std::array<float, 3> p0{};
        inputTriangulation_->getCellIncenter(cell.id_, cell.dim_, p0.data());

        points.push_back(p0[0]);
        points.push_back(p0[1]);
        points.push_back(p0[2]);

        points_pairOrigins.push_back(0);

        std::array<float, 3> p1{};
        inputTriangulation_->getCellIncenter(
          pairedCell.id_, pairedCell.dim_, p1.data());

        points.push_back(p1[0]);
        points.push_back(p1[1]);
        points.push_back(p1[2]);

        points_pairOrigins.push_back(1);

        cells.push_back(2);
        cells.push_back(pointId);
        cells.push_back(pointId + 1);

        cells_pairTypes.push_back(i);

        pointId += 2;
        cellId += 1;
      }
    }
  }

  numberOfPoints = pointId;
  numberOfCells = cellId;

  return 0;
}
