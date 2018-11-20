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

template<typename labelType> int finalize(
    vector<vector<Nodes>>& levelTimeNodesMap,
    vector<vector<Edges>>& levelTimeEdgesTMap,
    vector<vector<Edges>>& timeLevelEdgesNMap,
    int labelTypeId,
    string labelFieldName,

    vtkDataObject* trackingGraphObject
){
    auto trackingGraph = vtkUnstructuredGrid::SafeDownCast( trackingGraphObject );

    size_t nL = levelTimeNodesMap.size();
    size_t nT = levelTimeNodesMap[0].size();

    // Add Points
    {
        size_t nNodes = 0;
        for(size_t t=0; t<nT; t++)
            for(size_t l=0; l<nL; l++)
                nNodes += levelTimeNodesMap[l][t].size();

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
                for(auto& node: levelTimeNodesMap[l][t]){
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
        // Build node index offset vector
        vector<size_t> timeLevelOffsetMap(nT*nL+1);
        {
            timeLevelOffsetMap[0] = 0;
            size_t q = 1;
            for(size_t t=0; t<nT; t++)
                for(size_t l=0 ; l<nL; l++){
                    timeLevelOffsetMap[q] = timeLevelOffsetMap[q-1] + levelTimeNodesMap[l][t].size();
                    q++;
                }
        }

        size_t nEdgesT = 0;
        if(nT>1)
            for(size_t t=0; t<nT-1; t++)
                for(size_t l=0 ; l<nL; l++)
                    nEdgesT += levelTimeEdgesTMap[l][t].size()/3;

        size_t nEdgesN = 0;
        if(nL>1)
            for(size_t l=0; l<nL-1; l++)
                for(size_t t=0; t<nT; t++)
                    nEdgesN += timeLevelEdgesNMap[t][l].size()/3;

        auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
        cells->SetNumberOfValues( 3*nEdgesT + 3*nEdgesN );
        auto cellIds = (vtkIdType*) cells->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedLongLongArray> overlap = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        overlap->SetNumberOfValues( nEdgesT + nEdgesN );
        overlap->SetName("Overlap");
        auto overlapData = (unsigned long long*) overlap->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedCharArray> type = vtkSmartPointer<vtkUnsignedCharArray>::New();
        type->SetNumberOfValues( nEdgesT + nEdgesN );
        type->SetName("Type");
        auto typeData = (unsigned char*) type->GetVoidPointer(0);

        size_t q0=0, q1=0;

        // Tracking graphs
        if(nT>1)
            for(size_t t=1; t<nT; t++){
                for(size_t l=0; l<nL; l++){
                    auto& edges = levelTimeEdgesTMap[l][t-1];
                    for(size_t i=0, j=edges.size(); i<j; ){
                        cellIds[q0++] = 2;
                        cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ (t-1)*nL+l ] + edges[i++]);
                        cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ (t  )*nL+l ] + edges[i++]);
                        typeData[q1] = 0;
                        overlapData[q1++] = edges[i++];
                    }
                }
            }

        // Nesting trees
        if(nL>1)
            for(size_t l=1; l<nL; l++){
                for(size_t t=0; t<nT; t++){
                    auto& edges = timeLevelEdgesNMap[t][l-1];
                    size_t temp = t*nL;
                    for(size_t i=0, j=edges.size(); i<j; ){
                        cellIds[q0++] = 2;
                        cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ temp+(l-1) ] + edges[i++]);
                        cellIds[q0++] = (vtkIdType) (timeLevelOffsetMap[ temp+(l  ) ] + edges[i++]);
                        typeData[q1] = 1;
                        overlapData[q1++] = edges[i++];
                    }
                }
            }

        auto cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->SetCells(nEdgesT + nEdgesN, cells);
        trackingGraph->SetCells(VTK_LINE, cellArray);

        trackingGraph->GetCellData()->AddArray( overlap );
        trackingGraph->GetCellData()->AddArray( type );
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
    auto levels = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    /* Enforce following vtkMultiBlockDataSet structure:
        {
            level_0: {
                time_0: vtkPointSet,
                ...
                time_nT: vtkPointSet
            },
            ...
            level_nL: {
                time_0: vtkPointSet,
                ...
                time_nT: vtkPointSet
            }
        }
    */
    {
        // Check inputObject depth: 2->(level->time), 1->(-,time), 0->(-,-)
        auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast( inputObject );
        auto inputAsPS = vtkPointSet::SafeDownCast( inputObject );
        bool error = false;
        if(inputAsMB){
            size_t n = inputAsMB->GetNumberOfBlocks();

            // Check if blocks are vtkPointSets or vtkMultiBlockDataSets ...
            size_t psCounter = 0;
            size_t mbCounter = 0;
            for(size_t i=0; i<n; i++){
                auto block = inputAsMB->GetBlock(i);
                auto blockAsMB = vtkMultiBlockDataSet::SafeDownCast( block );
                auto blockAsPS = vtkPointSet::SafeDownCast( block );
                if(blockAsMB) mbCounter++;
                if(blockAsPS) psCounter++;
            }

            if(mbCounter==n){ // input already contains timesteps per level -> nothing to do
                levels->ShallowCopy( inputAsMB );
            } else if (psCounter==n){ // if input is a single list of vtkPointSets over time
                auto level = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                for(size_t i=0; i<n; i++){
                    level->SetBlock(
                        i,
                        vtkPointSet::SafeDownCast( inputAsMB->GetBlock(i) )
                    );
                }
                levels->SetBlock(0, level);
            } else { // Unexpected input structure
                error = true;
            }
        } else if(inputAsPS) {
            auto level = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            auto time = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            time->SetBlock(0, inputAsPS);
            level->SetBlock(0, time);
            levels->SetBlock(0, level);
        } else { // Unexpected input structure
            error = true;
        }

        if(error){
            dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Unable to convert input into 'vtkPointSet' collection.\n", fatalMsg);
            return 0;
        }
    }

    // Check input integrity
    this->LabelDataType = -1;
    size_t nT = 0;
    size_t nL = 0;
    {
        nL = levels->GetNumberOfBlocks();

        for(size_t l=0; l<nL; l++){
            auto times = vtkMultiBlockDataSet::SafeDownCast( levels->GetBlock(l) );
            size_t n = times->GetNumberOfBlocks();
            if(nT==0) nT = n;
            if(nT!=n) {
                dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Timeseries have unequal length.\n", fatalMsg);
                return 0;
            }

            for(size_t t=0; t<nT; t++){
                auto pointSet = vtkPointSet::SafeDownCast( times->GetBlock(t) );
                size_t nPoints = pointSet->GetNumberOfPoints();
                auto labels = pointSet->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );

                if(nPoints>0 && labels==nullptr){
                    dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' not found.\n" , fatalMsg);
                    return 0;
                }
                if(labels==nullptr) continue;

                int labelDataType = labels->GetDataType();
                if(this->LabelDataType<0) this->LabelDataType = labelDataType;
                if(this->LabelDataType != labelDataType){
                    dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels do not have same type across point sets.\n" , fatalMsg);
                    return 0;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Process Input
    // -------------------------------------------------------------------------

    Timer t;
    double t0 = 0;

    // Containers for nodes and edges
    vector<vector<Nodes>> levelTimeNodesMap(nL); // N
    vector<vector<Edges>> levelTimeEdgesTMap(nL); // E_T
    vector<vector<Edges>> timeLevelEdgesNMap(nT); // E_N

    // Reusable variables
    vtkPointSet* pointSet0=nullptr;
    vtkPointSet* pointSet1=nullptr;
    vtkAbstractArray* labels0=nullptr;
    vtkAbstractArray* labels1=nullptr;

    // Lambda function to fetch a specific block
    auto getData = [](vtkMultiBlockDataSet* mb, size_t time, size_t level, string labelFieldName, vtkPointSet*& pointSet, vtkAbstractArray*& labels){
        auto timesteps = vtkMultiBlockDataSet::SafeDownCast( mb->GetBlock( level ) );
        pointSet = vtkPointSet::SafeDownCast( timesteps->GetBlock( time ) );
        labels = pointSet->GetPointData()->GetAbstractArray( labelFieldName.data() );
    };

    // Compute Nodes
    {
        dMsg(cout, "[ttkTrackingFromOverlap] =======================================================\n" , infoMsg);
        dMsg(cout, "[ttkTrackingFromOverlap] Computing nodes\n" , infoMsg);
        t0 = t.getElapsedTime();

        for(size_t l=0; l<nL; l++){
            {
                stringstream msg;
                msg << "[ttkTrackingFromOverlap] -------------------------------------------------------" << endl
                    << "[ttkTrackingFromOverlap] Level Index: " << l << endl;
                dMsg(cout, msg.str(), infoMsg);
            }

            vector<Nodes>& timeNodesMap = levelTimeNodesMap[ l ];
            timeNodesMap.resize(nT);

            for(size_t t=0; t<nT; t++){
                getData(levels, t, l, this->GetLabelFieldName(), pointSet0, labels0);

                size_t nPoints0 = pointSet0->GetNumberOfPoints();
                if(nPoints0<1) continue;

                switch( this->LabelDataType ){
                    vtkTemplateMacro({
                        this->trackingFromOverlap.identifyNodes<VTK_TT>(
                            (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
                            (VTK_TT*) labels0->GetVoidPointer(0),
                            nPoints0,
                            timeNodesMap[ t ]
                        );
                    });
                }
            }
        }

        {
            stringstream msg;
            msg << "[ttkTrackingFromOverlap] ======================================================="<<endl
                << "[ttkTrackingFromOverlap] Nodes computed in " << (t.getElapsedTime()-t0) << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
        }

    }

    // Compute Tracking Graphs
    if(nT>1){
        dMsg(cout, "[ttkTrackingFromOverlap] =======================================================\n" , infoMsg);
        dMsg(cout, "[ttkTrackingFromOverlap] Computing tracking graphs\n" , infoMsg);
        t0 = t.getElapsedTime();

        for(size_t l=0; l<nL; l++){
            {
                stringstream msg;
                msg << "[ttkTrackingFromOverlap] -------------------------------------------------------" << endl
                    << "[ttkTrackingFromOverlap] Level Index: " << l << endl;
                dMsg(cout, msg.str(), infoMsg);
            }

            vector<Edges>& timeEdgesTMap = levelTimeEdgesTMap[ l ];
            timeEdgesTMap.resize(nT-1);

            for(size_t t=1; t<nT; t++){
                getData(levels, t-1, l, this->GetLabelFieldName(), pointSet0, labels0);
                getData(levels, t  , l, this->GetLabelFieldName(), pointSet1, labels1);

                size_t nPoints0 = pointSet0->GetNumberOfPoints();
                size_t nPoints1 = pointSet1->GetNumberOfPoints();
                if(nPoints0<1 || nPoints1<1) continue;

                switch( this->LabelDataType ){
                    vtkTemplateMacro({
                        this->trackingFromOverlap.computeOverlap<VTK_TT>(
                            (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
                            (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
                            (VTK_TT*) labels0->GetVoidPointer(0),
                            (VTK_TT*) labels1->GetVoidPointer(0),
                            nPoints0,
                            nPoints1,

                            timeEdgesTMap[ t-1 ]
                        );
                    });
                }
            }
        }

        {
            stringstream msg;
            msg << "[ttkTrackingFromOverlap] ======================================================="<<endl
                << "[ttkTrackingFromOverlap] Tracking graphs computed in " << (t.getElapsedTime()-t0) << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    // Compute Nesting Trees
    if(nL>1){
        dMsg(cout, "[ttkTrackingFromOverlap] =======================================================\n" , infoMsg);
        dMsg(cout, "[ttkTrackingFromOverlap] Computing nesting trees\n" , infoMsg);
        t0 = t.getElapsedTime();

        for(size_t t=0; t<nT; t++){
            {
                stringstream msg;
                msg << "[ttkTrackingFromOverlap] -------------------------------------------------------" << endl
                    << "[ttkTrackingFromOverlap] Time Index: " << t << endl;
                dMsg(cout, msg.str(), infoMsg);
            }

            vector<Edges>& levelEdgesNMap = timeLevelEdgesNMap[ t ];
            levelEdgesNMap.resize(nL-1);

            for(size_t l=1; l<nL; l++){
                getData(levels, t, l-1, this->GetLabelFieldName(), pointSet0, labels0);
                getData(levels, t, l  , this->GetLabelFieldName(), pointSet1, labels1);

                size_t nPoints0 = pointSet0->GetNumberOfPoints();
                size_t nPoints1 = pointSet1->GetNumberOfPoints();
                if(nPoints0<1 || nPoints1<1) continue;

                switch( this->LabelDataType ){
                    vtkTemplateMacro({
                        this->trackingFromOverlap.computeOverlap<VTK_TT>(
                            (float*)  pointSet0->GetPoints()->GetVoidPointer(0),
                            (float*)  pointSet1->GetPoints()->GetVoidPointer(0),
                            (VTK_TT*) labels0->GetVoidPointer(0),
                            (VTK_TT*) labels1->GetVoidPointer(0),
                            nPoints0,
                            nPoints1,

                            levelEdgesNMap[ l-1 ]
                        );
                    });
                }
            }
        }

        {
            stringstream msg;
            msg << "[ttkTrackingFromOverlap] ======================================================="<<endl
                << "[ttkTrackingFromOverlap] Nesting trees computed in " << (t.getElapsedTime()-t0) << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    // Get Output
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto trackingGraph = outInfo->Get(vtkDataObject::DATA_OBJECT());


    // Create Tracking Graph
    {
        dMsg(cout, "[ttkTrackingFromOverlap] =======================================================\n" , infoMsg);
        dMsg(cout, "[ttkTrackingFromOverlap] Meshing nested tracking graph\n" , infoMsg);
        t0 = t.getElapsedTime();

        switch( this->LabelDataType ){
            vtkTemplateMacro({
                finalize<VTK_TT>(
                    levelTimeNodesMap,
                    levelTimeEdgesTMap,
                    timeLevelEdgesNMap,
                    this->LabelDataType,
                    this->GetLabelFieldName(),
                    trackingGraph
                );
            });
        }

        {
            stringstream msg;
            msg << "[ttkTrackingFromOverlap] ======================================================="<<endl
                << "[ttkTrackingFromOverlap] Nested tracking graph meshed in " << (t.getElapsedTime()-t0) << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }
    {
        stringstream msg;
        msg << "[ttkTrackingFromOverlap] ======================================================="<<endl
            << "[ttkTrackingFromOverlap] Nested tracking graph generated in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 1;
}