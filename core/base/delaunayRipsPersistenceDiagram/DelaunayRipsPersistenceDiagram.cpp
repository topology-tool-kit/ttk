#include <DelaunayRipsPersistenceDiagram.h>

ttk::DelaunayRipsPersistenceDiagram::DelaunayRipsPersistenceDiagram() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("DelaunayRipsPD");
}

int ttk::DelaunayRipsPersistenceDiagram::execute(
  const PointCloud &points,
  MultidimensionalDiagram &ph) const {

#ifdef TTK_ENABLE_CGAL
  const int dim = points[0].size();
  if (dim > TTK_DELAUNAY_MAXIMUM_DIMENSION) {
    printErr("Input dimension too large: " + std::to_string(dim) + ">" + std::to_string(TTK_DELAUNAY_MAXIMUM_DIMENSION));
    return 1;
  }
  else {
    if (dim == 2) {
      FastRipsPersistenceDiagram2 FRPD(points);
      FRPD.setDebugLevel(debugLevel_);
      FRPD.computeDelaunayRips0And1Persistence(ph);
    }
    else if (dim == 3) {
      gph::runDelaunayRipsPersistenceDiagram3(points, ph);
      ph[0].emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair
    }
    else {
      gph::tryDimensions(points, ph);
      ph[0].emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair
      for (auto &diag : ph) {
        for (auto &[b,d] : diag) {
          b.first = {-1};
          d.first = {-1};
        }
      }
    }
  }
#else
  TTK_FORCE_USE(points);
  TTK_FORCE_USE(ph);
  printErr("TTK was not compiled with CGAL:");
  printErr("this filter is not available.");
#endif

  return 0;
}

int ttk::DelaunayRipsPersistenceDiagram::execute(
  const PointCloud &points,
  MultidimensionalDiagram &ph,
  std::vector<Generator1> &generators1,
  std::vector<Generator2> &generators2) const {

#ifdef TTK_ENABLE_CGAL
  const int dim = points[0].size();
  if (dim > 3) {
    printErr("Input dimension too large: " + std::to_string(dim) + ">3");
    return 1;
  }
  else {
    if (dim == 2) {
      FastRipsPersistenceDiagram2 FRPD(points);
      FRPD.setDebugLevel(debugLevel_);
      FRPD.computeDelaunayRips0And1Persistence(ph);
      FRPD.exportRips1Generators(generators1);
    }
    else if (dim == 3) {
      gph::runDelaunayRipsPersistenceDiagram3(points, ph, generators1, generators2);
      ph[0].emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair
    }
  }
#else
  TTK_FORCE_USE(points);
  TTK_FORCE_USE(ph);
  TTK_FORCE_USE(generators1);
  TTK_FORCE_USE(generators2);
  printErr("TTK was not compiled with CGAL:");
  printErr("this filter is not available.");
#endif

  return 0;
}