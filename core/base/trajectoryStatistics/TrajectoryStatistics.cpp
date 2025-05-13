#include <TrajectoryStatistics.h>
#include <Triangulation.h>

ttk::TrajectoryStatistics::TrajectoryStatistics() {
  this->setDebugMsgPrefix("TrajectoryStatistics");
}

int ttk::TrajectoryStatistics::findSurface(
  ttk::SimplexId                          startId,
  std::vector<ttk::SimplexId>            &surfVertex,
  const std::vector<std::vector<double>> &vertexScalars,
  std::vector<char>                      &visited,
  const double                            threshold,
  int                                    frame,
  const ttk::AbstractTriangulation       *triangulation
) {
  surfVertex.clear();
  // pile de DFS
  std::vector<ttk::SimplexId> stack;
  stack.reserve(1024);
  stack.push_back(startId);

  const double maxVal = threshold * 1.2;
  //this->printMsg("critical value : " + std::to_string(threshold) + " new thresh : " + std::to_string(maxVal)); 
  bool anyAdded = false;

  while(!stack.empty()) {
    auto vId = stack.back(); stack.pop_back();
    //this->printMsg("nouveau tour = " + std::to_string(vId));
    if(visited[vId]) continue;
    visited[vId] = 1;

    double val = vertexScalars[vId][frame];
    if(val > maxVal) continue;

    // on accepte ce sommet
    surfVertex.push_back(vId);
    //this->printMsg("current size = " + std::to_string(surfVertex.size()));
    anyAdded = true;

    // empiler ses voisins non visités
    const int nNbrs = triangulation->getVertexNeighborNumber(vId);
    for(int j = 0; j < nNbrs; ++j) {
      ttk::SimplexId nbr{-1};
      triangulation->getVertexNeighbor(vId, j, nbr);
      if(!visited[nbr]) {
        stack.push_back(nbr);
      }
    }
  }

  return anyAdded ? 1 : 0;
}


int ttk::TrajectoryStatistics::execute(
                std::vector<std::vector<int>> &trajTime,         
                std::vector<std::vector<double>> &trajX,          
                std::vector<std::vector<double>> &trajY,         
                std::vector<std::vector<double>> &trajZ,          
                std::vector<std::vector<int>>    &trajVertexId, 
                std::vector<std::vector<double>> vertexScalars,
                std::vector<int> &startFrames,          
                std::vector<int> &endFrames,            
                std::vector<int> &durations,
                std::vector<double> &VX,
                std::vector<double> &VY,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMoy,
                ttk::AbstractTriangulation *triangulation)  {


    const int numTraj = static_cast<int>(trajTime.size());

    // StartFrame - End - Duration 
    
    for (int i=0; i < numTraj ; i++){
        startFrames[i] = trajTime[i][0];
        endFrames[i] = trajTime[i].back();
        durations[i] = endFrames[i] - startFrames[i];
    }

    // VX / VY en pixel/frame 

    const int z_translation = trajZ[0][1] - trajZ[0][0]; 

    for (int i=0; i < numTraj; i++){
        double sommeVX = 0.0;
        double sommeVY = 0.0;
        for (int j=1; j<static_cast<int>(trajX[i].size()); j++){

            double dx = trajX[i][j] - trajX[i][j-1];
            double dy = trajY[i][j] - trajY[i][j-1];
            double dt = (trajZ[i][j] - trajZ[i][j-1])/z_translation;
            
            if (dt != 0){
                sommeVX += dx/dt;
                sommeVY += dy/dt;
            }

        }
        double numPoints = static_cast<double>(trajX[i].size()) - 1.0;
        VX[i] = sommeVX/numPoints; // checkez valeur etranges ici
        VY[i] = sommeVY/numPoints;
    }
    
    // Surface
   

    std::vector<std::vector<std::vector<ttk::SimplexId>>> allVertexDebris(numTraj);
//    #ifdef TTK_ENABLE_OPENMP
//    #pragma omp parallel for num_threads(this->threadNumber_)
//    #endif
    for(int i = 0; i < static_cast<int>(numTraj); ++i) {
        const int trajSize = static_cast<int>(trajVertexId[i].size());
        std::vector<int> trajSurfaces(trajSize);

        for(int j = 0; j < trajSize; ++j) {
            const int frame    = trajTime[i][j];
            const ttk::SimplexId vid = static_cast<ttk::SimplexId>(trajVertexId[i][j]);

            const double thresh = vertexScalars[vid][frame];

            std::vector<char> visited(vertexScalars.size(), 0);
            std::vector<ttk::SimplexId> surfVertex;
            surfVertex.reserve(64); // hypothèse d'une taille moyenne

            findSurface(vid, surfVertex, vertexScalars, visited, thresh,frame,  triangulation);
            
            //this->printMsg("nombre de sommet trouvé = " + std::to_string(surfVertex.size()));
            //allVertexDebris[i][frame] = std::move(surfVertex);     
            trajSurfaces.push_back(surfVertex.size());
            
        }
        std::sort(trajSurfaces.begin(), trajSurfaces.end(), [&](int a, int b){
                return a < b;
        });
        this->printMsg("Traj : " + std::to_string(i) + " surf min = " + std::to_string(trajSurfaces[0]) + " surf max = " + std::to_string(trajSurfaces.back()));
    }

    this->printMsg("End base");
    
    return 1;

}

