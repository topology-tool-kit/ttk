#include <ttkTrackingFromOverlap.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTrackingFromOverlap)

// int process(vtkDataObject* dataObject0, vtkDataObject* dataObject1){
//     auto pointSet1 = vtkPointSet::SafeDownCast( inMB->GetBlock(i) );
//     if(pointSet1==nullptr){
//         dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted to 'vtkPointSet'.\n", fatalMsg);
//         return 0;
//     }

//     auto labels1 = pointSet1->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
//     if(labels0==nullptr){
//         dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' not found.\n" , fatalMsg);
//         return 0;
//     }

//     this->trackingFromOverlap.identifyNodes(
//         (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
//         (VTK_TT*) labels1->GetVoidPointer(0),
//         pointSet1->GetNumberOfPoints(),
//         timeNodesMap[ i ]
//     );

//     this->trackingFromOverlap.sortCoordinates(
//         (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
//         pointSet1->GetNumberOfPoints(),
//         sortedIndicies[i%2]
//     );

//     this->trackingFromOverlap.track<VTK_TT>(
//         (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
//         (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
//         (VTK_TT*) labels0->GetVoidPointer(0),
//         (VTK_TT*) labels1->GetVoidPointer(0),
//         sortedIndicies[(i+1)%2],
//         sortedIndicies[(i  )%2],
//         pointSet0->GetNumberOfPoints(),
//         pointSet1->GetNumberOfPoints(),

//         timeNodesMap[ i-1 ],
//         timeNodesMap[ i   ],
//         timeEdgesMap[ i-1 ]
//     );
// }

template<typename labelType> int finalize(
    vector<vector<Nodes>>& timeLevelNodesMap,
    vector<vector<Edges>>& timeLevelEdgesMap,
    vtkDataObject* trackingGraphObject,
    int labelTypeId,
    string labelFieldName
){
    auto trackingGraph = vtkUnstructuredGrid::SafeDownCast( trackingGraphObject );

    size_t nT = timeLevelNodesMap.size();
    size_t nL = timeLevelNodesMap[0].size();

    // Add Points
    {
        size_t nNodes = 0;
        for(size_t t=0; t<nT; t++)
            for(size_t l=0; l<nL; l++)
                nNodes += timeLevelNodesMap[t][l].size();

        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints( nNodes );
        auto pointCoords = (float*) points->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedIntArray> time = vtkSmartPointer<vtkUnsignedIntArray>::New();
        time->SetName("TimeIndex");
        time->SetNumberOfComponents(1);
        time->SetNumberOfValues( nNodes );
        auto timeData = (unsigned int*) time->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedIntArray> level = vtkSmartPointer<vtkUnsignedIntArray>::New();
        level->SetName("LevelIndex");
        level->SetNumberOfComponents(1);
        level->SetNumberOfValues( nNodes );
        auto levelData = (unsigned int*) level->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedLongLongArray> size = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        size->SetName("Size");
        size->SetNumberOfComponents(1);
        size->SetNumberOfValues( nNodes );
        auto sizeData = (unsigned long long*) size->GetVoidPointer(0);

        vtkSmartPointer<vtkDataArray> label = vtkSmartPointer<vtkDataArray>::Take(
            vtkDataArray::CreateDataArray( labelTypeId )
        );
        label->SetName( labelFieldName.data() );
        label->SetNumberOfComponents(1);
        label->SetNumberOfValues( nNodes );
        auto labelData = (labelType*) label->GetVoidPointer(0);

        size_t q1=0, q2=0;
        for(size_t t=0; t<nT; t++){
            for(size_t l=0 ; l<nL; l++){
                for(auto& node: timeLevelNodesMap[t][l]){
                    pointCoords[q1++] = node.x;
                    pointCoords[q1++] = node.y;
                    pointCoords[q1++] = node.z;

                    timeData[q2]  = (unsigned int)t;
                    levelData[q2] = (unsigned int)l;
                    sizeData[q2]  = node.size;
                    labelData[q2] = boost::get<labelType>( node.label );
                    q2++;
                }
            }
        }

        trackingGraph->SetPoints(points);

        auto pointData = trackingGraph->GetPointData();
        pointData->AddArray( time );
        pointData->AddArray( level );
        pointData->AddArray( size );
        pointData->AddArray( label );
    }

    // Add Cells
    {
        size_t nEdges = 0;
        vector<size_t> timeLevelOffsetMap(nT*nL+1);
        {
            timeLevelOffsetMap[0] = 0;
            size_t q = 1;
            for(size_t t=0; t<nT; t++)
                for(size_t l=0 ; l<nL; l++){
                    timeLevelOffsetMap[q] = timeLevelOffsetMap[q-1] + timeLevelNodesMap[t][l].size();
                    q++;
                }
            for(size_t t=0; t<nT-1; t++)
                for(size_t l=0 ; l<nL; l++)
                    nEdges += timeLevelEdgesMap[t][l].size()/3;
        }

        auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
        cells->SetNumberOfValues( 3*nEdges );
        auto cellIds = (vtkIdType*) cells->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedLongLongArray> overlap = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        overlap->SetNumberOfValues( nEdges );
        overlap->SetName("Overlap");
        auto overlapData = (unsigned long long*) overlap->GetVoidPointer(0);

        for(size_t t=0, q0=0, q1=0; t<nT-1; t++){
            for(size_t l=0; l<nL; l++){
                auto& edges = timeLevelEdgesMap[t][l];
                for(size_t i=0, j=edges.size(); i<j; ){
                    cellIds[q0++] = 2;
                    cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ (t  )*nL+l ] + edges[i++]);
                    cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ (t+1)*nL+l ] + edges[i++]);
                    overlapData[q1++] = edges[i++];
                }
            }
        }

        auto cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->SetCells(nEdges, cells);
        trackingGraph->SetCells(VTK_LINE, cellArray);

        trackingGraph->GetCellData()->AddArray( overlap );
    }

    return 1;
}

int ttkTrackingFromOverlap::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkTrackingFromOverlap] RequestData" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Set Wrapper
    this->trackingFromOverlap.setWrapper(this);

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto inputObject = inInfo->Get(vtkDataObject::DATA_OBJECT());

    // -------------------------------------------------------------------------
    // Prepare Input
    // -------------------------------------------------------------------------
    /* Enforce following vtkMultiBlockDataSet structure:
        {
            time0: {
                level0: vtkPointSet,
                level1: vtkPointSet,
                ...
                levelLn: vtkPointSet
            },
            ...
            timeTn: {
                level0: vtkPointSet,
                level1: vtkPointSet,
                ...
                levelLn: vtkPointSet
            }
        }
    */
    auto timesteps = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    this->LabelDataType = -1;
    size_t nT = 0;
    size_t nL = 0;
    {
        auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast( inputObject );

        if(inputAsMB){ // Input contains timesteps
            nT = inputAsMB->GetNumberOfBlocks();

            // Check if blocks are already levels...
            size_t mbCounter = 0;
            size_t psCounter = 0;
            for(size_t i=0; i<nT; i++){
                auto block = inputAsMB->GetBlock(i);
                auto blockAsMB = vtkMultiBlockDataSet::SafeDownCast( block );
                auto blockAsPS = vtkPointSet::SafeDownCast( block );
                if(blockAsMB){
                    mbCounter++;
                    nL = blockAsMB->GetNumberOfBlocks();
                }
                if(blockAsPS){

                    auto labels = blockAsPS->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
                    if(labels)
                        this->LabelDataType = labels->GetDataType();

                    psCounter++;
                }
            }

            if(mbCounter==nT){ // input already contains levels -> nothing to do
                timesteps->ShallowCopy( inputAsMB );
            } else if (psCounter==nT){ // if input is list of vtkPointSets
                nL = 1;
                for(size_t i=0; i<nT; i++){
                    auto levels = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                    auto blockAsPS = vtkPointSet::SafeDownCast( inputAsMB->GetBlock(i) );

                    auto labels = blockAsPS->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
                    if(labels)
                        this->LabelDataType = labels->GetDataType();

                    levels->SetBlock(0, blockAsPS);
                    timesteps->SetBlock(i, levels);
                }
            } else {
                dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted into 'vtkPointSet' list.\n", fatalMsg);
                return 0;
            }
        } else {
            nT = 1;
            nL = 1;

            auto levels = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            auto inputAsPS = vtkPointSet::SafeDownCast( inputObject );
            if(!inputAsPS){
                dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted into 'vtkPointSet' list.\n", fatalMsg);
                return 0;
            }

            auto labels = inputAsPS->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
            if(labels)
                this->LabelDataType = labels->GetDataType();

            levels->SetBlock(0, inputAsPS);
            timesteps->SetBlock(0, levels);
        }
    }

    cout<< nT << " " << nL << " " << this->LabelDataType << endl;

    // -------------------------------------------------------------------------
    // Process Input
    // -------------------------------------------------------------------------

    // Containers for nodes and edges
    vector<vector<Nodes>> timeLevelNodesMap(nT);
    vector<vector<Edges>> timeLevelEdgesMap(nT-1);

    switch( this->LabelDataType ){
        vtkTemplateMacro({

            // Compute Nodes
            for(size_t t=0; t<nT; t++){
                auto levels = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock(t) );

                vector<Nodes>& levelNodesMap = timeLevelNodesMap[ t ];
                levelNodesMap.resize(nL);

                for(size_t l=0; l<nL; l++){
                    auto pointSet = vtkPointSet::SafeDownCast( levels->GetBlock(l) );
                    auto labels = pointSet->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );

                    this->trackingFromOverlap.identifyNodes<VTK_TT>(
                        (float*)  pointSet->GetPoints()->GetVoidPointer(0),
                        (VTK_TT*) labels->GetVoidPointer(0),
                        pointSet->GetNumberOfPoints(),
                        levelNodesMap[ l ]
                    );
                }
            }

            // Compute Edges
            for(size_t t=1; t<nT; t++){
                auto levels0 = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock( t-1 ) );
                auto levels1 = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock( t   ) );

                vector<Edges>& levelEdgesMap = timeLevelEdgesMap[ t-1 ];
                levelEdgesMap.resize(nL);

                for(size_t l=0; l<nL; l++){
                    auto pointSet0 = vtkPointSet::SafeDownCast( levels0->GetBlock(l) );
                    auto pointSet1 = vtkPointSet::SafeDownCast( levels1->GetBlock(l) );
                    auto labels0 = pointSet0->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
                    auto labels1 = pointSet1->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );

                    this->trackingFromOverlap.computeOverlap<VTK_TT>(
                        (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
                        (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
                        (VTK_TT*) labels0->GetVoidPointer(0),
                        (VTK_TT*) labels1->GetVoidPointer(0),
                        pointSet0->GetNumberOfPoints(),
                        pointSet1->GetNumberOfPoints(),

                        levelEdgesMap[ l ]
                    );
                }
            }

            // Get Output
            vtkInformation* outInfo = outputVector->GetInformationObject(0);
            auto trackingGraph = outInfo->Get(vtkDataObject::DATA_OBJECT());

            // Create Tracking Graph
            finalize<VTK_TT>(
                timeLevelNodesMap,
                timeLevelEdgesMap,
                trackingGraph,
                this->LabelDataType,
                this->GetLabelFieldName()
            );

        });
    }

    //     vector<Nodes>& levelNodesMap = timeLevelNodesMap[ t ];
    //     vector<Edges>& levelEdgesMap = timeLevelEdgesMap[ t ];
    //     levelNodesMap.resize(nL);
    //     // levelEdgesMap.resize(nL-1);

    //     // Init First Timestep across Levels
    //     // if(t==0){
    //     //     auto levels0 = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock(0) );

    //     //     for(size_t l=0; l<nL; l++){
    //     //         auto pointSet0 = vtkMultiBlockDataSet::SafeDownCast( levels0->GetBlock(l) );

    //     //         auto labels1 = pointSet1->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );

    //     //         if(l==0){
    //     //             this->trackingFromOverlap.identifyNodes(
    //     //                 (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
    //     //                 (VTK_TT*) labels0->GetVoidPointer(0),
    //     //                 pointSet0->GetNumberOfPoints(),
    //     //                 timeNodesMap[ 0 ]
    //     //             );
    //     //         }
    //     //     }

    //     //     // auto pointSet0 = vtkMultiBlockDataSet::SafeDownCast( levels0->GetBlock(0) );
    //     // }



    //     auto levels0 = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock(t  ) );
    //     auto levels1 = vtkMultiBlockDataSet::SafeDownCast( timesteps->GetBlock(t+1) );

    //     for(size_t l=0; l<nL-1; l++){

    //         auto pointSet0 = vtkMultiBlockDataSet::SafeDownCast( levels0->GetBlock(l) );
    //         auto pointSet1 = vtkMultiBlockDataSet::SafeDownCast( levels1->GetBlock(l) );

    //         auto labels0 = pointSet0->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
    //         auto labels1 = pointSet1->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );

    //         if(l==0){
    //             this->trackingFromOverlap.identifyNodes(
    //                 (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
    //                 (VTK_TT*) labels0->GetVoidPointer(0),
    //                 pointSet0->GetNumberOfPoints(),
    //                 timeNodesMap[ 0 ]
    //             );
    //         }
    //     }
    // }

    // double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    // double nIterations = inMB->GetInformation()->Get( vtkDataObject::DATA_TIME_STEP() );

    // bool iterativeMode = nIterations!=0;
    // bool iterativeMode = false;

    // if(!iterativeMode && nBlocks>0){
    //     vtkSmartPointer<vtkMultiBlockDataSet> inMB

    //     // Containers for nodes and edges
    //     vector<Edges> timeEdgesMap(nBlocks-1);
    //     vector<Nodes> timeNodesMap(nBlocks);

    //     // Containers for sorted indicies
    //     vector<vector<size_t>> sortedIndicies(2);

    //     this->LabelDataType = -1;

    //     auto pointSet0 = vtkPointSet::SafeDownCast( inMB->GetBlock(0) );
    //     if(pointSet0==nullptr){
    //         dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted to 'vtkPointSet'.\n", fatalMsg);
    //         return 0;
    //     }

    //     auto labels0 = pointSet0->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
    //     if(labels0==nullptr){
    //         dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' not found.\n" , fatalMsg);
    //         return 0;
    //     }

    //     switch( labels0->GetDataType() ){
    //         ttkTemplateMacro({

    //             // Identify nodes of first point set
    //             this->trackingFromOverlap.identifyNodes(
    //                 (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
    //                 (VTK_TT*) labels0->GetVoidPointer(0),
    //                 pointSet0->GetNumberOfPoints(),
    //                 timeNodesMap[ 0 ]
    //             );

    //             // Sort coordinates of first point set
    //             this->trackingFromOverlap.sortCoordinates(
    //                 (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
    //                 pointSet0->GetNumberOfPoints(),
    //                 sortedIndicies[0]
    //             );

    //             // Iterate over remaining blocks
    //             for(size_t i=1; i<nBlocks; i++){
    //                 auto pointSet1 = vtkPointSet::SafeDownCast( inMB->GetBlock(i) );
    //                 if(pointSet1==nullptr){
    //                     dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted to 'vtkPointSet'.\n", fatalMsg);
    //                     return 0;
    //                 }

    //                 auto labels1 = pointSet1->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
    //                 if(labels0==nullptr){
    //                     dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' not found.\n" , fatalMsg);
    //                     return 0;
    //                 }

    //                 this->trackingFromOverlap.identifyNodes(
    //                     (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
    //                     (VTK_TT*) labels1->GetVoidPointer(0),
    //                     pointSet1->GetNumberOfPoints(),
    //                     timeNodesMap[ i ]
    //                 );

    //                 this->trackingFromOverlap.sortCoordinates(
    //                     (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
    //                     pointSet1->GetNumberOfPoints(),
    //                     sortedIndicies[i%2]
    //                 );

    //                 this->trackingFromOverlap.track<VTK_TT>(
    //                     (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
    //                     (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
    //                     (VTK_TT*) labels0->GetVoidPointer(0),
    //                     (VTK_TT*) labels1->GetVoidPointer(0),
    //                     sortedIndicies[(i+1)%2],
    //                     sortedIndicies[(i  )%2],
    //                     pointSet0->GetNumberOfPoints(),
    //                     pointSet1->GetNumberOfPoints(),

    //                     timeNodesMap[ i-1 ],
    //                     timeNodesMap[ i   ],
    //                     timeEdgesMap[ i-1 ]
    //                 );

    //                 pointSet0 = pointSet1;
    //                 labels0 = labels1;

    //                 this->updateProgress( ((float)i)/((float)(nBlocks)) );
    //             }

    //             // Get Output
    //             vtkInformation* outInfo = outputVector->GetInformationObject(0);
    //             auto trackingGraph = outInfo->Get(vtkDataObject::DATA_OBJECT());

    //             // Create Tracking Graph
    //             finalize<VTK_TT>(
    //                 timeNodesMap,
    //                 timeEdgesMap,
    //                 trackingGraph,
    //                 labels0->GetDataType(),
    //                 this->GetLabelFieldName()
    //             );

    //             // Print Status
    //             dMsg(cout, "[ttkTrackingFromOverlap] Tracking Graph Finalized\n", infoMsg);
    //             this->updateProgress( 1 );
    //         });
    //     }
    // }

    return 1;
}