#include <PlanarGraphLayout.h>

#if TTK_ENABLE_GRAPHVIZ
#include <graphviz/cgraph.h>
#include <graphviz/gvc.h>
#endif

// Compute Dot Layout
int ttk::PlanarGraphLayout::computeDotLayout(
  // Input
  const vector<size_t> &nodeIndicies,
  const string &dotString,

  // Output
  float *layout) const {
#if TTK_ENABLE_GRAPHVIZ
  Timer t;

  dMsg(cout, "[ttkPlanarGraphLayout] Computing layout ... ", timeMsg);

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
    Agnode_t *n = agnode(G, const_cast<char *>(to_string(i).data()), 0);
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

  // ---------------------------------------------------------
  // Print status
  // ---------------------------------------------------------
  {
    stringstream msg;
    msg << "done (" << t.getElapsedTime() << " s)" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
#else
  dMsg(cout,
       "[ttkPlanarGraphLayout] ERROR: This filter requires GraphViz to compute "
       "a layout.\n",
       fatalMsg);
  return 0;
#endif
}
