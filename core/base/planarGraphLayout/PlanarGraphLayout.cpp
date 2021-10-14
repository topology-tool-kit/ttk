#include <PlanarGraphLayout.h>

#if TTK_ENABLE_GRAPHVIZ
#include <cgraph.h>
#include <gvc.h>
#endif

ttk::PlanarGraphLayout::PlanarGraphLayout() {
  this->setDebugMsgPrefix("PlanarGraphLayout");
}
ttk::PlanarGraphLayout::~PlanarGraphLayout() {
}

// Compute Dot Layout
int ttk::PlanarGraphLayout::computeDotLayout(
  // Output
  float *layout,

  // Input
  const std::vector<size_t> &nodeIndicies,
  const std::string &dotString) const {
#ifdef TTK_ENABLE_GRAPHVIZ
  Timer t;

  this->printMsg("Computing layout", 0, debug::LineMode::REPLACE);

  // ---------------------------------------------------------
  // Init GraphViz
  // ---------------------------------------------------------
  Agraph_t *G = agmemread(dotString.data());
  GVC_t *gvc = gvContext();
  gvLayout(gvc, G, "dot");

  // ---------------------------------------------------------
  // Get layout data from GraphViz
  // ---------------------------------------------------------
  for(auto i : nodeIndicies) {
    Agnode_t *n = agnode(G, const_cast<char *>(std::to_string(i).data()), 0);
    if(n != nullptr) {
      auto &coord = ND_coord(n);
      size_t offset = i * 2;
      layout[offset] = coord.x / 72; // points to inches
      layout[offset + 1] = coord.y / 72; // points to inches
    }
  }

  // ---------------------------------------------------------
  // Free GraphViz memory
  // ---------------------------------------------------------
  gvFreeLayout(gvc, G);
  agclose(G);
  gvFreeContext(gvc);

  this->printMsg("Computing layout", 1, t.getElapsedTime());

  return 1;
#else
  TTK_FORCE_USE(layout);
  TTK_FORCE_USE(nodeIndicies);
  TTK_FORCE_USE(dotString);

  this->printErr("This filter requires GraphViz to compute a layout.");
  return 0;
#endif // TTK_ENABLE_GRAPHVIZ
}
