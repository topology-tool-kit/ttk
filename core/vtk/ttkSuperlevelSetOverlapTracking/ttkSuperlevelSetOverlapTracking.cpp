#include                  <ttkSuperlevelSetOverlapTracking.h>

#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkStringArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkCellData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkSuperlevelSetOverlapTracking)

int ttkSuperlevelSetOverlapTracking::doIt(vtkTable* inputTable, vtkUnstructuredGrid* outputGrid){

    Memory m;

    superlevelSetOverlapTracking_.setWrapper(this);

    auto indexPathMap = vtkStringArray::SafeDownCast(inputTable->GetColumnByName( "path" ));
    auto indexTimeMap = vtkIntArray::SafeDownCast(inputTable->GetColumnByName( "time" ));

    size_t tn = indexTimeMap->GetNumberOfValues();

    vtkSmartPointer<vtkXMLImageDataReader> vtkReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    vtkImageData* image0data;
    vtkImageData* image1data;

    vtkReader->SetFileName( indexPathMap->GetValue(0).c_str() );
    vtkReader->Update();
    image0data = vtkReader->GetOutput();

    int* dim = image0data->GetDimensions();
    superlevelSetOverlapTracking_.setDim(dim);

    vector< vector<Node> > timeNodeMap(tn);
    vector< vector<Edge> > timeEdgeMap(tn-1);

    auto createLabelArray = [](int* dim){
        vtkSmartPointer<vtkFloatArray> labels = vtkSmartPointer<vtkFloatArray>::New();
        labels->SetNumberOfComponents(1);
        labels->SetName("Label");
        labels->SetNumberOfValues( dim[0]*dim[1]*dim[2] );
        return labels;
    };

    auto writeImageData = [](string path, vtkImageData* imageData, vtkSmartPointer<vtkFloatArray> labels){
        return;
        cout<<"Writing: "<< path <<endl;

        vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        vtkWriter->SetFileName(path.c_str());

        imageData->GetPointData()->AddArray( labels );
        vtkWriter->SetInputData(imageData);
        vtkWriter->Write();
    };

    vtkSmartPointer<vtkFloatArray> image0Labels = createLabelArray(dim);
    vtkSmartPointer<vtkFloatArray> image1Labels = createLabelArray(dim);

    switch(image0data->GetPointData()->GetArray(0)->GetDataType()){
        vtkTemplateMacro({

            auto computeLabels = [](string path, VTK_TT level, void* data, void* labels, vector<Node>& nodes, float time, SuperlevelSetOverlapTracking& superlevelSetOverlapTracking){
                cout<<"Labeling "<< path << ": " << flush;
                superlevelSetOverlapTracking.computeLabels<VTK_TT>(
                    level,
                    (VTK_TT*) data,
                    (float*) labels,
                    nodes,
                    time
                );
                cout<< "#" << nodes.size() << endl;
            };

            computeLabels(
                indexPathMap->GetValue(0),
                (VTK_TT) this->Level,
                image0data->GetPointData()->GetArray(0)->GetVoidPointer(0),
                image0Labels->GetVoidPointer(0),
                timeNodeMap[0],
                0,
                superlevelSetOverlapTracking_
            );
            writeImageData( indexPathMap->GetValue(0) + ".labels.vti", image0data, image0Labels);

            for(size_t t=0; t<tn-1; t++){
                bool flip = t%2==0;
                auto l0 = flip ? image0Labels : image1Labels;
                auto l1 = flip ? image1Labels : image0Labels;

                vtkReader->SetFileName( indexPathMap->GetValue(t+1).c_str() );
                vtkReader->Update();
                image1data = vtkReader->GetOutput();

                computeLabels(
                    indexPathMap->GetValue(t+1),
                    (VTK_TT) this->Level,
                    image1data->GetPointData()->GetArray(0)->GetVoidPointer(0),
                    l1->GetVoidPointer(0),
                    timeNodeMap[t+1],
                    t+1,
                    superlevelSetOverlapTracking_
                );
                writeImageData( indexPathMap->GetValue(t+1) + ".labels.vti", image1data, l1);

                cout<<"Tracking\n\t"<< indexPathMap->GetValue(t) << "\n\t" << indexPathMap->GetValue(t+1) <<endl;
                superlevelSetOverlapTracking_.execute<VTK_TT>(
                    (float*) l0->GetVoidPointer(0),
                    (float*) l1->GetVoidPointer(0),
                    timeEdgeMap[t]
                );
                cout<<"Done"<<endl;
            }
        });
    }

    // cout<< "Computing layout..." << flush;
    // superlevelSetOverlapTracking_.computeLayout(timeNodeMap, timeEdgeMap);
    // cout<< "done" << endl;

    cout<< "Computing Branches..." << flush;
    superlevelSetOverlapTracking_.computeBranches(timeNodeMap, timeEdgeMap);
    cout<< "done" << endl;

    cout<< "Generating Grid..." << flush;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();


    // Point Data
    {
        vtkSmartPointer<vtkDoubleArray> size = vtkSmartPointer<vtkDoubleArray>::New();
        size->SetNumberOfComponents(1);
        size->SetName("Size");

        vtkSmartPointer<vtkUnsignedLongLongArray> branches = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        branches->SetNumberOfComponents(1);
        branches->SetName("NodeBranch");

        size_t q=0;
        for(size_t t=0; t<tn; t++){
            auto& nodes = timeNodeMap[t];
            for(size_t i=0; i<nodes.size(); i++){
                auto& n = nodes[i];
                points->InsertNextPoint( n.x, n.y, n.z );
                size->InsertTuple1(q, (double) n.size);
                branches->InsertTuple1(q++, n.branch);
            }
        }

        mesh->SetPoints(points);
        mesh->GetPointData()->AddArray(size);
        mesh->GetPointData()->AddArray(branches);
    }

    // Cell Data
    {
        vtkSmartPointer<vtkUnsignedLongLongArray> branches = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        branches->SetNumberOfComponents(1);
        branches->SetName("EdgeBranch");

        vector<size_t> offset(tn+1);
        offset[0] = 0;
        for(size_t t=1; t<tn; t++){
            offset[t] = offset[t-1] + (&timeNodeMap[t-1])->size();
        }

        size_t q=0;
        for(size_t t=0; t<tn-1; t++){
            auto& edges = timeEdgeMap[t];
            for(auto& e: edges){
                vtkIdType ids[2] = { (vtkIdType)(e.n0+offset[t]), (vtkIdType)(e.n1+offset[t+1]) };
                mesh->InsertNextCell(VTK_LINE, 2, ids);
                branches->InsertTuple1(q++, e.branch);
            }
        }
        mesh->GetCellData()->AddArray(branches);
    }

    // Field Data
    {
        // auto paths = vtkSmartPointer<vtkStringArray>::New();
        // paths->SetNumberOfComponents(1);
        // paths->SetName("Files");
        // paths->Resize(tn);

        auto times = vtkSmartPointer<vtkDoubleArray>::New();
        times->SetNumberOfComponents(1);
        times->SetName("Time");
        times->Resize(tn);
        for(size_t t=0; t<tn; t++){
            times->SetValue(t, (double) indexTimeMap->GetValue(t));
        }
        mesh->GetFieldData()->AddArray( times );
    }

    outputGrid->ShallowCopy(mesh);
    cout<< "done" << endl;

    {
        stringstream msg;
        msg << "[ttkSuperlevelSetOverlapTracking] Memory usage: " << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
}

int ttkSuperlevelSetOverlapTracking::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid* outputGrid = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  doIt(
      inputTable,
      outputGrid
  );

  return 1;
}