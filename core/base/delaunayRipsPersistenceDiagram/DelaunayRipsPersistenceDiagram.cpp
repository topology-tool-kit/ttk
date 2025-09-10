#include <DelaunayRipsPersistenceDiagram.h>

ttk::DelaunayRipsPersistenceDiagram::DelaunayRipsPersistenceDiagram() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("DelaunayRipsPD");
}

int ttk::DelaunayRipsPersistenceDiagram::execute(
  const rpd::PointCloud &points,
  rpd::MultidimensionalDiagram &ph,
  std::vector<rpd::Generator> &generators) const {

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
      if(OutputGenerators)
        FRPD.exportRips1Generators(generators);
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

  for (auto &diag : ph) {
    for (auto &[b,d] : diag)
      std::cout << b.second << " " << d.second << std::endl;
    std::cout << "----------" << std::endl;
  }
#else
  printErr("TTK was not compiled with CGAL:");
  printErr("this filter is not available.");
#endif

  return 0;
}