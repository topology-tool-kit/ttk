#include "Graph.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include<iostream>
#endif

using namespace std;
using namespace ttk;
using namespace ftr;

Graph::Graph()
{
}

Graph::~Graph()
{
}

void Graph::alloc()
{
#ifndef TTK_ENABLE_KAMIKAZE
   if (nbVerts_ == nullVertex) {
      cout << "[FTR Graph]: ERROR, setNumberOfVertices not called before alloc in Graph" << endl;
   }
#endif
   leaves_.reserve(nbVerts_/3);
   nodes_.reserve(nbVerts_/2);
   arcs_.reserve(nbVerts_/2);
   segmentation_.resize(nbVerts_);
}

void Graph::init()
{
   fillVector<graphElmt>(segmentation_, nullGraphElmt);
}
