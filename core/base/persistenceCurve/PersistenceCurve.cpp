#include <Geometry.h>
#include <PersistenceCurve.h>

ttk::PersistenceCurve::PersistenceCurve() {
  setDebugMsgPrefix("PersistenceCurve");
}

int ttk::PersistenceCurve::execute(std::array<PlotType, 4> &plots,
                                   const DiagramType &diagram) const {

  const auto epsilon{Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS)};

  // copy input diagram and sort the copy by pair persistence
  auto sortedDiagram{diagram};

  TTK_PSORT(this->threadNumber_, sortedDiagram.begin(), sortedDiagram.end(),
            [](const PersistencePair &a, const PersistencePair &b) {
              return a.persistence() < b.persistence();
            });

  for(const auto &pair : sortedDiagram) {
    // per dimension
    plots[pair.dim].emplace_back(
      std::max(pair.persistence(), epsilon), plots[pair.dim].size());
    // all pairs
    plots[3].emplace_back(
      std::max(pair.persistence(), epsilon), plots[3].size());
  }

  for(auto &plot : plots) {
    for(auto &el : plot) {
      el.second = plot.size() - el.second;
    }
  }

  // look at the first finite pair of dimension 1
  const auto firstPairDim1 = std::find_if(
    diagram.begin(), diagram.end(),
    [](const PersistencePair &p) { return p.dim == 1 && p.isFinite; });

  // if critical type of death is maximum, the diagram was computed on
  // a 2D dataset
  const auto datasetIs2D
    = firstPairDim1 != diagram.end()
      && firstPairDim1->death.type == ttk::CriticalType::Local_maximum;

  // enforce plots[2] is for saddle-max pairs
  if(datasetIs2D) {
    plots[2] = std::move(plots[1]);
  }

  return 0;
}
