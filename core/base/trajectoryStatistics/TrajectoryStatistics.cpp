#include <TrajectoryStatistics.h>
#include <Triangulation.h>

ttk::TrajectoryStatistics::TrajectoryStatistics() {
  this->setDebugMsgPrefix("TrajectoryStatistics");
}


int ttk::TrajectoryStatistics::execute(
                std::vector<std::vector<int>> &trajTime,         
                std::vector<std::vector<double>> &trajX,          
                std::vector<std::vector<double>> &trajY,         
                std::vector<std::vector<double>> &trajZ,          
                std::vector<std::vector<int>>    &trajVertexId,  
                std::vector<int> &startFrames,          
                std::vector<int> &endFrames,            
                std::vector<int> &durations,
                std::vector<double> &VX,
                std::vector<double> &VY,
                ttk::AbstractTriangulation *triangulation)  {

    // StartFrame - End - Duration 
    
    for (int i=0; i<trajTime.size(); i++){
        startFrames[i] = trajTime[i][0];
        endFrames[i] = trajTime[i].back();
        durations[i] = endFrames[i] - startFrames[i];
    }

    // VX / VY en pixel/frame 

    const int z_translation = trajZ[0][1] - trajZ[0][0]; 

    for (int i=0; i<trajX.size(); i++){
        double sommeVX, sommeVY = 0.0;
        for (int j=1; j<trajX[i].size(); j++){
            double dx = trajX[i][j] - trajX[i][j-1];
            double dy = trajY[i][j] - trajY[i][j-1];
            double dt = (trajZ[i][j] - trajZ[i][j-1])/z_translation;
            
            if (dt != 0){
                sommeVX += dx/dt;
                sommeVY += dy/dt;
            }
        }
        VX[i] = sommeVX/trajX[i].size();
        VY[i] = sommeVY/trajY[i].size();
    }
    

    
    return 1;

}
