/// \ingroup base
/// \class ttk::TopoMap
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023
///
///
/// This module defines the %TopoMap class that embeds points into 2D. The input
/// is a distance matrix representing the distances between any pair of points
/// in high dimension. It computes the 2D coordinates of a projection of the
/// points which preserves the distances of a minimum spanning tree computed by
/// Kruskal's algorithm. It also preserves the components Kruskal's algorithm
/// builds iteratively.
///
///
/// \b Related \b publication: \n
/// "Topomap: A 0-dimensional homology preserving projection of high-dimensional
/// data"\n Harish Doraiswamy, Julien Tierny, Paulo J. S. Silva, Luis Gustavo
/// Nonato, and Claudio Silva\n Proc. of IEEE VIS 2020.\n IEEE Transactions on
/// Visualization and Computer Graphics 27(2): 561-571, 2020.

/// \sa TopoMap

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <UnionFind.h>

// STL includes
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <utility> // For std::pair
#include <vector>

// cos and sin computations, so we are cautious on testing equalities between
// double.
static double Epsilon{ttk::Geometry::powInt(10.0, -FLT_DIG + 2)};
static double EpsilonDBL{ttk::Geometry::powInt(10.0, -DBL_DIG + 2)};

// Normalizes a given vector.
template <typename T>
static void
  computeUnitVector(const T *coordOrig, const T *coordDest, T *coordVect);

// Rotates the first argument by angle around center.
template <typename T>
static void rotate(T *ptToRotate, const T *center, double angle);

// Rotate the set of points by the given angle, considering the given center as
// center of the rotation.
template <typename T>
static void
  rotatePolygon(std::vector<T> &coords, T *centerCoords, const double angle);

// Real call to the convex hull engine: either Qhull or Boost.
bool computeConvexHull_aux(const std::vector<double> &coords,
                           std::vector<size_t> &res,
                           std::string &errMsg);

namespace ttk {

  /**
   * The TopoMap class computes a 2D projection of the points given as input,
   * either by their distance matrix or by their coordinates.
   */
  class TopoMap : virtual public Debug {

  public:
    /** We rely on constructing a minimal spanning tree. We allow for the two
     * main algorithms to build one: Kruskal and Prim.
     */
    enum class STRATEGY { KRUSKAL = 0, PRIM = 1 };

    TopoMap();
    TopoMap(size_t angularSampleNb, bool checkMST, STRATEGY strategy)
      : AngularSampleNb(angularSampleNb), CheckMST(checkMST),
        Strategy(strategy) {

      // inherited from Debug: prefix will be printed at the beginning of every
      // msg
      this->setDebugMsgPrefix("TopoMap");
    }
    ~TopoMap() override;

    /**
     * @brief Computes the TopoMap projection
     *
     * @param[out] outputCoords the final coordinates of the points, computed by
     TopoMap
     *
     * @param[out] insertionTime an optional array indicating for each point at
     what time it was inserted in the projection.
     *
     * @param[in] inputMatrix the high-dimension input distance matrix or
     coordinates of the points
     *
     * @param[in] isDistMat equal to true if inputMatrix is a distance matrix,
     false if it stores input point coordinates
     *
     * @param[in] n the number of input points

     * @return 0 in case of success.
     */
    template <typename T>
    int execute(T *outputCoords,
                int *insertionTime,
                const std::vector<T> &inputMatrix,
                bool isDistMat,
                size_t n);

  protected:
    size_t AngularSampleNb{2};
    bool CheckMST{false};
    bool errorConvexHull{false};
    STRATEGY Strategy{STRATEGY::KRUSKAL};

  private:
    // Fills its last two arguments with the position of the end of the previous
    // (resp. next) edge so that idCenter is the point at the angle at the
    // intersection of the previous and next edges.
    template <typename T>
    bool getPrevNextEdges(const std::vector<size_t> &idsPtsPolygon,
                          size_t idCenter,
                          const T *allCoords,
                          T *coordPrev,
                          T *coordPost) const;

    // Tries to find the best angle of rotation for the two components. Updates
    // the coordiates of their vertices accordingly.
    template <typename T>
    T rotateMergingCompsBest(const std::vector<size_t> &hull1,
                             const std::vector<size_t> &hull2,
                             const std::vector<size_t> &comp1,
                             const std::vector<size_t> &comp2,
                             size_t iPt1,
                             size_t iPt2,
                             const std::vector<T> &distMatrix,
                             T *allCoords,
                             size_t n,
                             size_t angularSampleNb);

    template <typename T>
    bool computeConvexHull(T *allCoords,
                           const std::vector<size_t> &compPtsIds,
                           std::vector<size_t> &idsInHull) const;
  };

  // Sketch of the algorithm:
  // 1. Sort the edges
  // 2. For each edge uv with cost c_uv
  //      a. Get connected components comp(u), comp(v)
  //
  //      b. Compute convex hull of comp(u) and of (comp(v))
  //
  //      c. Find in hull(u) (resp. hull(v)) the point pU (resp. pV) in the
  //      convex hull which is closest to u (resp. v) in high dimension, u if
  //      possible.
  //
  //      d. Translate comp(u) such that pU lies on the line of the bissector of
  //      angleV, at distance c_uv of pV, and on the opposite direction of
  //      hull(v). angleU (resp. angleV) is the angle of hull(u) (resp. hull(v))
  //      at vertex pU (resp. pV).
  //
  //      e. Rotate comp(u) such that the bissector of angleU prolongates the
  //      one of angleU, and such that pU and pV are closest than any pair of
  //      points in comp(u),comp(v): the two components are facing each other at
  //      pU and pV.
  //
  //
  //      f. Try several rotations of comp(u) and comp(v) such that pU and pV
  //      still are the closest pair of vertices from the two components. We
  //      keep the angle of rotation which minimizes the difference between the
  //      high dimension distance matrix and the new distance matrix between the
  //      two components.
  //
  //      g. comp(u) and comp(v) were merged, we now consider them as a single
  //      component.

  template <typename T>
  int ttk::TopoMap::execute(T *outputCoords,
                            int *insertionTime,
                            const std::vector<T> &inputMatrix,
                            bool isDistMat,
                            size_t n) {
    ttk::Timer timer;

#ifdef TTK_ENABLE_QHULL
#ifndef Qhull_FOUND
    this->printWrn(
      "Qhull was enabled but it is not installed or was not found.");
    this->printWrn("Defaulting to Boost support instead.");
#endif
#endif

#ifndef Qhull_FOUND
    this->printMsg("Using Boost for convex hulls.");
    this->printWrn("Bugs have been reported in Boost's implementation.");
    this->printWrn("Consider enabling Qhull instead.");
#endif

    if(AngularSampleNb < 2) {
      this->printWrn("The number of angular samples set is less than 2.");
      this->printWrn("Setting it to 2 (minimum possible value).");
      this->AngularSampleNb = 2;
    }
    // We first set the minimum value we consider as null. We multiply by 100
    // the smallest significative value because we use sqrt/cos/sin so we take
    // precautionary measures.
    if(std::is_same<T, double>::value) {
      Epsilon = ttk::Geometry::powInt(10.0, -DBL_DIG + 2);
    } else if(std::is_same<T, float>::value) {
      Epsilon = ttk::Geometry::powInt(10.0, -FLT_DIG + 2);
    } else {
      printErr("Error, template type must be either double or float.\n");
      return 1;
    }

    if(isDistMat) {
      this->printMsg("Input data: " + std::to_string(n) + " points.");
    } else {
      this->printMsg("Input data: " + std::to_string(n) + " points ("
                     + std::to_string(inputMatrix.size() / n)
                     + " dimensions).");
    }

    std::vector<T> computedDistMatrix(n * n);
    if(!isDistMat) {
      size_t dim = inputMatrix.size() / n;
      if(n * dim != inputMatrix.size()) {
        this->printErr("Error, the coordinates input matrix has invalid size.");
        return 1;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) \
  shared(computedDistMatrix)
#endif
      for(size_t i1 = 0; i1 < n; i1++) {
        for(size_t i2 = i1; i2 < n; i2++) {
          T dist = ttk::Geometry::distance<T>(
            &inputMatrix[dim * i1], &inputMatrix[dim * i2], dim);
          computedDistMatrix[i1 * n + i2] = dist;
          computedDistMatrix[i2 * n + i1] = dist;
        }
      }
    } else {
      if(inputMatrix.size() != n * n) {
        this->printErr("Invalid size for the distance matrix.");
        return 1;
      }
    }

    std::vector<double> edgesMSTBefore,
      edgesMSTAfter; // Only used for minimum spanning tree checking;
    edgesMSTBefore.reserve(n);
    const std::vector<T> &distMatrix
      = isDistMat ? inputMatrix : computedDistMatrix;

    if(insertionTime != nullptr)
      for(size_t i = 0; i < n; i++)
        insertionTime[i] = -1;

    // 1. Sorting the edges
    using heapType = std::pair<T, std::pair<size_t, size_t>>;
    std::priority_queue<heapType, std::vector<heapType>, std::greater<heapType>>
      edgeHeap;
    if(this->Strategy == STRATEGY::KRUSKAL) {
      for(size_t u1 = 0; u1 < n; u1++) {
        for(size_t u2 = u1 + 1; u2 < n; u2++) {
          edgeHeap.push(
            std::make_pair(distMatrix[u1 * n + u2], std::make_pair(u1, u2)));
        }
      }
    } else if(this->Strategy == STRATEGY::PRIM) {
      for(size_t u1 = 0; u1 < n; u1++) {
        edgeHeap.push(std::make_pair(distMatrix[u1], std::make_pair(0, u1)));
      }
    } else {
      this->printErr("Invalid stategy for the MST: only Kruskal (0) or Prim "
                     "(1) algorithm can be used.");
      return 1;
    }

    // ufVectors is the union find vector
    std::map<UnionFind *, std::vector<size_t>> ufToComp;
    std::vector<UnionFind> ufVector(n);
    for(size_t i = 0; i < n; i++) {
      ufToComp[&ufVector[i]].emplace_back(i);
    }

    for(size_t i = 0; i < 2 * n; i++) {
      outputCoords[i] = 0;
    }

    T finalDistortion = 0;
    // 2. Processing the edges
    size_t nbEdgesMerged = 0;
    double MSTWeight = 0;

    std::set<size_t>
      stillToDo; // List of vertices not yet embedded, used for Prim strategy.
    if(this->Strategy == STRATEGY::PRIM) {
      for(size_t i = 1; i < n; i++) {
        stillToDo.insert(i);
      }
    }
    while(nbEdgesMerged != n - 1) {
      if(edgeHeap.empty()) {
        this->printErr("Error building the MST. Aborting.");
        return 1;
      }
      const auto &elt = edgeHeap.top();
      T edgeCost = elt.first;
      size_t u = elt.second.first;
      size_t v = elt.second.second;
      edgeHeap.pop();

      // 2.a Getting the components
      UnionFind *reprU = ufVector[u].find();
      UnionFind *reprV = ufVector[v].find();
      if(reprU == reprV) { // Already in the same component
        continue;
      }
      nbEdgesMerged++;
      MSTWeight += edgeCost;

      if(insertionTime != nullptr) {
        if(insertionTime[u] == -1) {
          insertionTime[u] = nbEdgesMerged;
        }
        if(insertionTime[v] == -1) {
          insertionTime[v] = nbEdgesMerged;
        }
      }
      edgesMSTBefore.push_back(edgeCost);
      std::vector<size_t> &compU = ufToComp[reprU], &compV = ufToComp[reprV];
      size_t idSmall = compU.size() < compV.size() ? u : v;
      size_t idBig = idSmall == u ? v : u;
      std::vector<size_t> &compSmall = idSmall == u ? compU : compV;
      std::vector<size_t> &compBig = idSmall == u ? compV : compU;
      size_t nBig = compBig.size();
      size_t nSmall = compSmall.size();
      std::vector<size_t> idsInHullSmall, idsInHullBig;

      // 2.b Computing the convex hull.
      // We retrieve the current coordinates of the big component..
      bool statusHull = true;
      statusHull
        = statusHull | computeConvexHull(outputCoords, compBig, idsInHullBig);
      statusHull = statusHull
                   | computeConvexHull(outputCoords, compSmall, idsInHullSmall);
      if(!statusHull) {
        this->printErr("Aborting.");
        return 1;
      }

      size_t idChosenBig = idsInHullBig[0], idChosenSmall = idsInHullSmall[0];

      // 2.c We want to select, among all vertices in the convex hull, the one
      // which is closest to the vertex of the edge we work on.
      for(size_t vert : idsInHullBig) {
        T dist = distMatrix[vert * n + idBig];
        if(dist < distMatrix[idChosenBig * n + idBig]) {
          idChosenBig = vert;
        }
      }

      for(size_t vert : idsInHullSmall) {
        T dist = distMatrix[vert * n + idSmall];
        if(dist < distMatrix[idChosenSmall * n + idSmall]) {
          idChosenSmall = vert;
        }
      }

      size_t sizeBigHull = idsInHullBig.size(),
             sizeSmallHull = idsInHullSmall.size();
      std::vector<T> coordsBigHull(sizeBigHull * 2),
        coordsSmallHull(sizeSmallHull * 2);
      for(size_t iHull = 0; iHull < sizeBigHull; iHull++) {
        size_t vert = idsInHullBig[iHull];
        coordsBigHull[iHull * 2] = outputCoords[vert * 2];
        coordsBigHull[iHull * 2 + 1] = outputCoords[vert * 2 + 1];
      }
      for(size_t iHull = 0; iHull < sizeSmallHull; iHull++) {
        size_t vert = idsInHullSmall[iHull];
        coordsSmallHull[iHull * 2] = outputCoords[vert * 2];
        coordsSmallHull[iHull * 2 + 1] = outputCoords[vert * 2 + 1];
      }

      // Identifying the angles we are working on from the convex hulls.
      T coordPrevBig[2], coordPostBig[2];
      T coordPrevSmall[2], coordPostSmall[2];
      if(!getPrevNextEdges(idsInHullSmall, idChosenSmall, outputCoords,
                           coordPrevSmall, coordPostSmall)) {
        return 1;
      }
      if(!getPrevNextEdges(idsInHullBig, idChosenBig, outputCoords,
                           coordPrevBig, coordPostBig)) {
        return 1;
      }
      T coordPtSmall[2] = {
        outputCoords[2 * idChosenSmall], outputCoords[2 * idChosenSmall + 1]};
      T coordPtBig[2]
        = {outputCoords[2 * idChosenBig], outputCoords[2 * idChosenBig + 1]};

      // Computing the angles.
      double angleSmall = ttk::Geometry::angle2DUndirected(
        coordPtSmall, coordPrevSmall, coordPostSmall);
      double angleBig = ttk::Geometry::angle2DUndirected(
        coordPtBig, coordPrevBig, coordPostBig);
      if(angleSmall - M_PI > EpsilonDBL || angleBig - M_PI > EpsilonDBL) {
        this->printErr("Error, angle out of bound (greater than pi).");
      }
      T coordscenterBig[2] = {coordPrevBig[0], coordPrevBig[1]};
      T coordscenterSmall[2] = {coordPrevSmall[0], coordPrevSmall[1]};
      // Computing the coordinates of the bisectors.
      rotate(coordscenterSmall, coordPtSmall, angleSmall / 2);
      rotate(coordscenterBig, coordPtBig, angleBig / 2);

      T unitcenterBigVect[2], unitcenterSmallVect[2];
      computeUnitVector(
        &outputCoords[idChosenBig * 2], coordscenterBig, unitcenterBigVect);
      computeUnitVector(&outputCoords[idChosenSmall * 2], coordscenterSmall,
                        unitcenterSmallVect);

      // Computing the new coordinates for the chosen point on the small
      // components, and then translating the whole component.
      T goalCoordChosenSmall[2]
        = {outputCoords[idChosenBig * 2] - edgeCost * unitcenterBigVect[0],
           outputCoords[idChosenBig * 2 + 1] - edgeCost * unitcenterBigVect[1]};
      if(sizeBigHull == 1) {
        goalCoordChosenSmall[0] = outputCoords[idChosenBig * 2] + edgeCost;
        goalCoordChosenSmall[1] = outputCoords[idChosenBig * 2 + 1];
      }
      T smallCompMoveVect[2]
        = {goalCoordChosenSmall[0] - outputCoords[idChosenSmall * 2],
           goalCoordChosenSmall[1] - outputCoords[idChosenSmall * 2 + 1]};

      T distBaryPointSmall = ttk::Geometry::distance2D(
        &outputCoords[idChosenSmall * 2], coordscenterSmall);
      T preFinalPosBarySmall[2] = {coordscenterSmall[0] + smallCompMoveVect[0],
                                   coordscenterSmall[1] + smallCompMoveVect[1]};
      T finalPosBarySmall[2]
        = {goalCoordChosenSmall[0] - unitcenterBigVect[0] * distBaryPointSmall,
           goalCoordChosenSmall[1] - unitcenterBigVect[1] * distBaryPointSmall};

      // Performing the translation + rotation such that the bisectors in
      // compSmall and compBig are aligned, and the two components are facing
      // each other, not crossing.
      for(size_t curIdSmall : compSmall) {
        outputCoords[curIdSmall * 2] += smallCompMoveVect[0];
        outputCoords[curIdSmall * 2 + 1] += smallCompMoveVect[1];

        // 2.e Performing the rotation.
        double rotationAngle = ttk::Geometry::angle2DUndirected(
          goalCoordChosenSmall, preFinalPosBarySmall, finalPosBarySmall);
        if(nSmall > 1 && std::isfinite(rotationAngle)) {
          rotate(
            &outputCoords[curIdSmall * 2], goalCoordChosenSmall, rotationAngle);
        }
      }

      // 2.f Trying several rotations of the two components and keeping the one
      // which makes the new distance matrix closest to the one provided in the
      // input.
      if(nBig > 1) {
        finalDistortion = rotateMergingCompsBest(
          idsInHullSmall, idsInHullBig, compSmall, compBig, idChosenSmall,
          idChosenBig, distMatrix, outputCoords, n, this->AngularSampleNb);
        if(finalDistortion < -1) {
          return 1;
        }
      }

      T finalDist = ttk::Geometry::distance2D(
        &outputCoords[2 * idChosenSmall], &outputCoords[2 * idChosenBig]);
      if(fabs(finalDist - edgeCost) > Epsilon) {
        this->printErr(
          "The distance we set is too far from the goal distance.");
      }

      UnionFind *unionRepr = UnionFind::makeUnion(reprU, reprV);
      UnionFind *otherRepr = (unionRepr == reprU) ? reprV : reprU;
      if(unionRepr != reprU && unionRepr != reprV) {
        printErr("Error with the union find module results. Aborting.");
        return 1;
      }
      std::vector<size_t> &unionSet = ufToComp[unionRepr];
      std::vector<size_t> &otherSet = ufToComp[otherRepr];

      unionSet.insert(unionSet.end(), otherSet.begin(), otherSet.end());

      ufToComp.erase(otherRepr);

      if(this->Strategy == STRATEGY::PRIM) {
        stillToDo.erase(v);
        for(size_t uToDo : stillToDo) {
          edgeHeap.push(std::make_pair(
            distMatrix[v * n + uToDo], std::make_pair(v, uToDo)));
        }
      }
    }
    std::stringstream ssDistortion;
    ssDistortion << std::scientific << 2 * finalDistortion;
    this->printMsg("Non normalized distance matrix distortion: "
                   + ssDistortion.str());
    this->printMsg(
      "(for normalized distortion values, use distanceMatrixDistorsion).");

    if(this->errorConvexHull) {
      this->printWrn(
        "Warning, somme convex hull computations were wrong, so the projection "
        "may be less optimal than it could have been. If you are using Boost, "
        "consider switching to Qhull to avoid this.");
    }
    this->printMsg("The minimum spanning tree has weight "
                   + std::to_string(MSTWeight) + ".");

    // We check that the lengths of the edges selected to build a minimum
    // spanning tree are preserved by our algorithm, as should be.
    if(CheckMST) {

      edgesMSTAfter.reserve(n);
      this->printMsg("Checking the new minimum spanning tree.");
      std::priority_queue<heapType, std::vector<heapType>,
                          std::greater<heapType>>
        edgeHeapAfter;
      if(this->Strategy == STRATEGY::KRUSKAL) {
        for(size_t u1 = 0; u1 < n; u1++) {
          for(size_t u2 = u1 + 1; u2 < n; u2++) {
            edgeHeapAfter.push(
              std::make_pair(ttk::Geometry::distance2D(
                               &outputCoords[2 * u1], &outputCoords[2 * u2]),
                             std::make_pair(u1, u2)));
          }
        }
      } else if(this->Strategy == STRATEGY::PRIM) {
        for(size_t u1 = 0; u1 < n; u1++) {
          edgeHeapAfter.push(std::make_pair(
            ttk::Geometry::distance2D(&outputCoords[0], &outputCoords[2 * u1]),
            std::make_pair(0, u1)));
        }
      }

      std::map<UnionFind *, std::vector<size_t>> ufToCompAfter;
      std::vector<UnionFind> ufVectorAfter(n);
      std::vector<UnionFind *> ufPtrVectorAfter(n);
      for(size_t i = 0; i < n; i++) {
        ufPtrVectorAfter[i] = &ufVectorAfter[i];
        ufToCompAfter[ufPtrVectorAfter[i]].push_back(i);
      }

      nbEdgesMerged = 0;
      stillToDo.clear();
      for(size_t i = 1; i < n; i++)
        stillToDo.insert(i);
      while(nbEdgesMerged != n - 1) {
        if(edgeHeapAfter.empty()) {
          this->printErr("Error building the MST for the check. Aborting.");
          return 1;
        }
        const auto &elt = edgeHeapAfter.top();
        T edgeCost = elt.first;
        size_t u = elt.second.first;
        size_t v = elt.second.second;
        edgeHeapAfter.pop();

        UnionFind *reprU = ufPtrVectorAfter[u]->find();
        UnionFind *reprV = ufPtrVectorAfter[v]->find();
        if(reprU == reprV) { // Already in the same component
          continue;
        }
        UnionFind *unionRepr = UnionFind::makeUnion(reprU, reprV);
        UnionFind *otherRepr = (unionRepr == reprU) ? reprV : reprU;
        if(unionRepr != reprU && unionRepr != reprV) {
          printErr("Error with the union find module results. Aborting.");
          return 1;
        }
        nbEdgesMerged++;
        std::vector<size_t> &unionSet = ufToCompAfter[unionRepr];
        std::vector<size_t> &otherSet = ufToCompAfter[otherRepr];

        unionSet.insert(unionSet.end(), otherSet.begin(), otherSet.end());

        ufPtrVectorAfter[u] = unionRepr;
        ufPtrVectorAfter[v] = unionRepr;

        ufToCompAfter.erase(otherRepr);
        if(this->Strategy == STRATEGY::PRIM) {
          stillToDo.erase(v);
          for(size_t uToDo : stillToDo) {
            edgeHeapAfter.push(
              std::make_pair(ttk::Geometry::distance2D(
                               &outputCoords[2 * v], &outputCoords[2 * uToDo]),
                             std::make_pair(v, uToDo)));
          }
        }

        edgesMSTAfter.push_back(edgeCost);
      }

      const T epsilonCheck = Epsilon * 5;
      double sumDiff = 0, maxDiff = 0;
      size_t nbProblematic = 0;
      if(this->Strategy == STRATEGY::PRIM) {
        sort(edgesMSTBefore.begin(), edgesMSTBefore.end());
        sort(edgesMSTAfter.begin(), edgesMSTAfter.end());
      }
      for(size_t i = 0; i < edgesMSTBefore.size(); i++) {
        if(fabs(edgesMSTBefore[i] - edgesMSTAfter[i]) >= epsilonCheck) {
          nbProblematic++;
          double diff = fabs(edgesMSTBefore[i] - edgesMSTAfter[i]);
          sumDiff += diff;
          maxDiff = std::max(diff, maxDiff);
        }
      }
      if(nbProblematic != 0) {
        this->printWrn("Error on " + std::to_string(nbProblematic)
                       + " edge(s) of the MST.");

        std::stringstream ssMax, ssAvg;
        ssMax << std::scientific << maxDiff;
        ssAvg << std::scientific << sumDiff / nbProblematic;
        this->printWrn("The maximum error is " + ssMax.str()
                       + " and the average (on the number of errors) error is "
                       + ssAvg.str() + ".");
      }
    }

    this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
    this->printMsg("Complete", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator

    return 0;
  }

  template <typename T>
  bool TopoMap::getPrevNextEdges(const std::vector<size_t> &idsPtsPolygon,
                                 size_t idCenter,
                                 const T *allCoords,
                                 T *coordPrev,
                                 T *coordPost) const {
    size_t n = idsPtsPolygon.size();
    size_t iPtPrev = 0, iPtPost = 0;
    bool found = false;

    for(size_t i = 0; i < n; i++) {
      if(idsPtsPolygon[i] == idCenter) {
        found = true;
        iPtPost = idsPtsPolygon[(i + 1) % n];
        iPtPrev = idsPtsPolygon[(i + n - 1) % n];
        break;
      }
    }
    if(!found) {
      printErr("Error, we could not find the edges incident to the point we "
               "chose for the rotation of the component. Aborting.");
      return false;
    }

    coordPrev[0] = allCoords[2 * iPtPrev];
    coordPrev[1] = allCoords[2 * iPtPrev + 1];
    coordPost[0] = allCoords[2 * iPtPost];
    coordPost[1] = allCoords[2 * iPtPost + 1];

    return true;
  }

  template <typename T>
  T TopoMap::rotateMergingCompsBest(const std::vector<size_t> &hull1,
                                    const std::vector<size_t> &hull2,
                                    const std::vector<size_t> &comp1,
                                    const std::vector<size_t> &comp2,
                                    size_t iPt1,
                                    size_t iPt2,
                                    const std::vector<T> &distMatrix,
                                    T *allCoords,
                                    size_t n,
                                    size_t angularSampleNb) {
    // The distance between the two components.
    T shortestDistPossible
      = ttk::Geometry::distance2D(&allCoords[2 * iPt1], &allCoords[2 * iPt2]);
    T coordPt1[2] = {allCoords[2 * iPt1], allCoords[2 * iPt1 + 1]};
    T coordPt2[2] = {allCoords[2 * iPt2], allCoords[2 * iPt2 + 1]};
    size_t comp1Size = comp1.size(), comp2Size = comp2.size();

    T coordPrev1[2], coordPost1[2];
    T coordPrev2[2], coordPost2[2];
    getPrevNextEdges(hull1, iPt1, allCoords, coordPrev1, coordPost1);
    getPrevNextEdges(hull2, iPt2, allCoords, coordPrev2, coordPost2);

    double angle1
      = ttk::Geometry::angle2DUndirected(coordPt1, coordPrev1, coordPost1);
    double angle2
      = ttk::Geometry::angle2DUndirected(coordPt2, coordPrev2, coordPost2);
    if(angle1 - M_PI > EpsilonDBL || angle2 - M_PI > EpsilonDBL) {
      this->printErr("One angle of the convex hull is greater than pi. "
                     "Convexity error, aborting.");
      return -2;
    }
    T coordBissect1[2] = {coordPrev1[0], coordPrev1[1]};
    T coordBissect2[2] = {coordPrev2[0], coordPrev2[1]};
    rotate(coordBissect1, coordPt1, angle1 / 2);
    rotate(coordBissect2, coordPt2, angle2 / 2);

    // If a component has only one vertex, no need to try several angles.
    double semiAngle1 = comp1Size == 1 ? M_PI / 2 : angle1 / 2;
    double semiAngle2 = comp2Size == 1 ? M_PI / 2 : angle2 / 2;

    double angleMax1 = M_PI / 2 - semiAngle1, angleMin1 = -angleMax1;
    double angleMax2 = M_PI / 2 - semiAngle2, angleMin2 = -angleMax2;
    double step1 = (angleMax1 - angleMin1) / angularSampleNb,
           step2 = (angleMax2 - angleMin2) / angularSampleNb;
    double bestAnglePair[2] = {0, 0};
    T bestScore = 3.e38; // Near the maximum value for float

    std::vector<std::vector<T>> origDistMatrix(comp1Size);
    for(size_t i = 0; i < comp1Size; i++) {
      origDistMatrix[i].resize(comp2Size);
      for(size_t j = 0; j < comp2Size; j++) {
        origDistMatrix[i][j] = distMatrix[comp1[i] * n + comp2[j]];
      }
    }
    std::vector<T> initialCoords1(2 * comp1Size), initialCoords2(2 * comp2Size);
    for(size_t i = 0; i < comp1.size(); i++) {
      size_t iPt = comp1[i];
      initialCoords1[2 * i] = allCoords[2 * iPt];
      initialCoords1[2 * i + 1] = allCoords[2 * iPt + 1];
    }
    for(size_t i = 0; i < comp2.size(); i++) {
      size_t iPt = comp2[i];
      initialCoords2[2 * i] = allCoords[2 * iPt];
      initialCoords2[2 * i + 1] = allCoords[2 * iPt + 1];
    }

    size_t nbIter1 = std::isfinite(step1) ? angularSampleNb + 1 : 1;
    size_t nbIter2 = std::isfinite(step2) ? angularSampleNb + 1 : 1;

    if(step1 * nbIter1 < 0.001) { // No need to split such a small angle
      nbIter1 = 1;
    }
    if(step2 * nbIter2 < 0.001) { // No need to split such a small angle
      nbIter2 = 1;
    }

    bool errorDuringLoop = false;
    bool foundValidAngle = false;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) shared(allCoords)
#endif
    // We enumerate the rotation angles for the first component.
    for(size_t i1 = 0; i1 < nbIter1; i1++) {
      std::vector<T> coords1Test(2 * comp1Size), coords2Test(2 * comp2Size);
      coords1Test = initialCoords1;
      double testAngle1 = angleMin1 + step1 * i1;
      rotatePolygon(coords1Test, coordPt1, testAngle1);

      // For each first rotation angle, we enumerate the rotations angles for
      // the second componenent.
      for(size_t i2 = 0; i2 < nbIter2; i2++) {
        bool curAngleError = false;
        coords2Test = initialCoords2;
        double testAngle2 = angleMin2 + i2 * step2;
        rotatePolygon(coords2Test, coordPt2, testAngle2);

        // We compute a distortion score if we rotate the two components by the
        // current angles.
        T curScore = 0;
        for(size_t i = 0; i < comp1Size; i++) {
          T coordARotate[2] = {coords1Test[2 * i], coords1Test[2 * i + 1]};
          for(size_t j = 0; j < comp2Size; j++) {
            T coordBRotate[2] = {coords2Test[2 * j], coords2Test[2 * j + 1]};
            T newDist = ttk::Geometry::distance2D(coordARotate, coordBRotate);
            curScore += (newDist - origDistMatrix[i][j])
                        * (newDist - origDistMatrix[i][j]);
            if(newDist + Epsilon < shortestDistPossible) {
              curAngleError = true;
              errorDuringLoop = true;
              break;
            }
          }
        }

        if(curAngleError) {
          continue;
        }

        // Now if the score is the best so far, we keep this pair angle as the
        // best one we tried so far.
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          // In case of errors in convex hull computation, we still have a valid
          // angle (to preserve the MST).
          foundValidAngle = true;
          // This pair of angle minimises the distortion and the choice is
          // deterministic with threads.
          if(curScore < bestScore
             || (curScore == bestScore
                 && (testAngle1 < bestAnglePair[0]
                     || (testAngle1 == bestAnglePair[0]
                         && testAngle2 < bestAnglePair[1])))) {
            bestScore = curScore;
            bestAnglePair[0] = std::isfinite(testAngle1) ? testAngle1 : 0;
            bestAnglePair[1] = std::isfinite(testAngle2) ? testAngle2 : 0;
          }
        }
      }
    }

    // No valid angle were found, because the convex hull is invalid;
    if(!foundValidAngle) {
      this->printErr(
        "No valid angle was found, due to errors in the convex hull "
        "computation. Please consider using Qhull instead of Boost. Aborting.");
      return -10;
    }

    // Applying the chosen rotations to the two components.
    for(size_t i1 : comp1) {
      rotate(&allCoords[2 * i1], coordPt1, bestAnglePair[0]);
    }
    for(size_t i2 : comp2) {
      rotate(&allCoords[2 * i2], coordPt2, bestAnglePair[1]);
    }

    if(errorDuringLoop) {
      this->errorConvexHull = true;
    }

    return bestScore;
  }

  template <typename T>
  bool TopoMap::computeConvexHull(T *allCoords,
                                  const std::vector<size_t> &compPtsIds,
                                  std::vector<size_t> &idsInHull) const {
    size_t nbPoint = compPtsIds.size();
    std::vector<double> compCoords(2 * nbPoint);
    for(size_t i = 0; i < nbPoint; i++) {
      size_t vertId = compPtsIds[i];
      compCoords[2 * i] = static_cast<double>(allCoords[2 * vertId]);
      compCoords[2 * i + 1] = static_cast<double>(allCoords[2 * vertId + 1]);
    }

    if(nbPoint <= 2) {
      idsInHull.push_back(compPtsIds[0]);
      if(nbPoint == 2) {
        double dist = ttk::Geometry::distance2D(&compCoords[0], &compCoords[2]);

        if(dist > Epsilon) {
          idsInHull.push_back(compPtsIds[1]);
        }
      }
      return true;
    }

    // Testing if all points are colinear. If so, we take extremal ones. We do
    // this to avoid that the convex hull algorithm does not detect it and gives
    // a wrong result.
    double dirVect[2] = {0, 0};
    bool areColinear = true;

    size_t idFirstDistinct = 0;
    while(idFirstDistinct < nbPoint && fabs(dirVect[0]) < Epsilon
          && fabs(dirVect[1]) < Epsilon) {
      idFirstDistinct++;
      dirVect[0] = compCoords[2 * idFirstDistinct] - compCoords[0];
      dirVect[1] = compCoords[2 * idFirstDistinct + 1] - compCoords[1];
      if(fabs(dirVect[0]) < Epsilon) {
        dirVect[0] = 0;
      }
      if(fabs(dirVect[1]) < Epsilon) {
        dirVect[1] = 0;
      }
    }

    // Minimum/maximum values for x and y coordinates.
    int idMins[2] = {compCoords[0] < compCoords[2] ? 0 : 1,
                     compCoords[1] < compCoords[3] ? 0 : 1};
    int idMaxs[2] = {compCoords[0] > compCoords[2] ? 0 : 1,
                     compCoords[1] > compCoords[3] ? 0 : 1};

    const double *pt0 = &compCoords[0];
    const double *ptDistinct = &compCoords[2 * idFirstDistinct];

    for(size_t iPt = idFirstDistinct + 1; iPt < nbPoint; iPt++) {
      double curVect[2] = {compCoords[2 * iPt] - compCoords[0],
                           compCoords[2 * iPt + 1] - compCoords[1]};
      if(fabs(curVect[0]) < Epsilon && fabs(curVect[1]) < Epsilon) {
        continue;
      }
      const double *ptCur = &compCoords[2 * iPt];
      if(!ttk::Geometry::isTriangleColinear2D(
           pt0, ptDistinct, ptCur, EpsilonDBL)) {
        areColinear = false;
        break;
      }

      // Updating the min and max values for x and y coordinates.
      if(compCoords[2 * iPt] < compCoords[2 * idMins[0]]) {
        idMins[0] = iPt;
      }
      if(compCoords[2 * iPt + 1] < compCoords[2 * idMins[1]]) {
        idMins[1] = iPt;
      }
      if(compCoords[2 * iPt] > compCoords[2 * idMaxs[0]]) {
        idMaxs[0] = iPt;
      }
      if(compCoords[2 * iPt + 1] > compCoords[2 * idMaxs[1]]) {
        idMaxs[1] = iPt;
      }
    }

    if(areColinear) {
      if(fabs(dirVect[0]) > Epsilon) {
        idsInHull.push_back(compPtsIds[idMins[0]]);
        idsInHull.push_back(compPtsIds[idMaxs[0]]);
      } else {
        idsInHull.push_back(compPtsIds[idMins[1]]);
        idsInHull.push_back(compPtsIds[idMaxs[1]]);
      }
      return true;
    }

    // Computing the convex hull using Boost or Qhull according to the build.
    std::string errMsg;
    bool status = computeConvexHull_aux(compCoords, idsInHull, errMsg);
    if(!status) {
      this->printErr(errMsg);
      return false;
    }

    // In its C++ API, Qhull does not return the points of the convex hull in
    // (anti-)clockwise order. Boost does not seem consistent with the
    // orientation of the hull. So we reorder the points by sorting them
    // according to the angle they form with the center and a fictitious point
    // at its right.
    double sumX = 0, sumY = 0;
    for(size_t u : idsInHull) {
      sumX += compCoords[2 * u];
      sumY += compCoords[2 * u + 1];
    }
    double bary[2] = {sumX / idsInHull.size(), sumY / idsInHull.size()};
    double baryRight[2] = {bary[0] + 2, bary[1]};
    std::vector<std::pair<double, size_t>> ptsToSort;
    ptsToSort.reserve(idsInHull.size());
    for(size_t u : idsInHull) {
      const double curPt[2] = {compCoords[2 * u], compCoords[2 * u + 1]};
      double curAngle
        = ttk::Geometry::angle2DUndirected(bary, baryRight, curPt);
      ptsToSort.emplace_back(std::make_pair(-curAngle, u));
    }

    // We sort the points according to the angle, keeping the indices
    // associated.
    sort(ptsToSort.begin(), ptsToSort.end());
    for(size_t i = 0; i < ptsToSort.size(); i++)
      idsInHull[i] = ptsToSort[i].second;

    // We obtained the indices of the points belonging to the hull. We replace
    // them by the global ids of the chosen points.
    for(size_t &x : idsInHull) {
      x = compPtsIds[x];
    }
    return true;
  }
} // namespace ttk

template <typename T>
static void
  computeUnitVector(const T *coordOrig, const T *coordDest, T *coordVect) {
  T tmp[2] = {coordDest[0] - coordOrig[0], coordDest[1] - coordOrig[1]};
  T dist = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);
  if(dist < Epsilon) {
    coordVect[0] = 1;
    coordVect[1] = 0;
  } else {
    coordVect[0] = tmp[0] / dist;
    coordVect[1] = tmp[1] / dist;
  }
}

template <typename T>
static void rotate(T *ptToRotate, const T *center, double angle) {
  const double &xCtr = center[0], &yCtr = center[1];
  T &xPt = ptToRotate[0], &yPt = ptToRotate[1];
  const double dx = xPt - xCtr, dy = yPt - yCtr;
  xPt = dx * cos(angle) - dy * sin(angle) + xCtr;
  yPt = dx * sin(angle) + dy * cos(angle) + yCtr;
}

template <typename T>
static void
  rotatePolygon(std::vector<T> &coords, T *centerCoords, const double angle) {
  double xCenter = centerCoords[0], yCenter = centerCoords[1];
  size_t nbPoint = coords.size() / 2;
  for(size_t iPt = 0; iPt < nbPoint; iPt++) {
    T &x = coords[iPt * 2], &y = coords[iPt * 2 + 1];
    double xNew, yNew;
    xNew = (x - xCenter) * cos(angle) - (y - yCenter) * sin(angle) + xCenter;
    yNew = (y - yCenter) * cos(angle) + (x - xCenter) * sin(angle) + yCenter;
    x = xNew;
    y = yNew;
  }
}
