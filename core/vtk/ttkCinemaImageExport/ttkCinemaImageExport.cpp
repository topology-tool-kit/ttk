#include                  <ttkCinemaImageExport.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkDataSet.h>
#include                  <vtkMultiBlockDataSet.h>
#include                  <vtkWindowToImageFilter.h>
#include                  <vtkDataSetMapper.h>
#include                  <vtkCompositePolyDataMapper.h>
#include                  <vtkCompositePolyDataMapper2.h>
#include                  <vtkActor.h>
#include                  <vtkRenderer.h>
#include                  <vtkRenderWindow.h>
#include                  <vtkCamera.h>
#include                  <vtkCompositePolyDataMapper2.h>
#include                  <vtkPolyDataMapper.h>
#include                  <vtkCompositeDataDisplayAttributes.h>

#include                  <vtkPointSet.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkPointData.h>
#include                  <vtkFieldData.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkMath.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaImageExport)

int ttkCinemaImageExport::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaImageExport] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;
    Timer t;
    double t0=0;

    vtkInformation* inputGeomertyInfo = inputVector[0]->GetInformationObject(0);
    auto inputGeomerty = inputGeomertyInfo->Get(vtkDataObject::DATA_OBJECT());

    vtkInformation* inGridInfo = inputVector[1]->GetInformationObject(0);
    auto inputGrid = vtkPointSet::SafeDownCast( inGridInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Visualize
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    if(inputGeomerty->IsA("vtkMultiBlockDataSet")){
        auto x = vtkMultiBlockDataSet::SafeDownCast(inputGeomerty)->GetBlock(0);
        mapper->SetInputData( vtkUnstructuredGrid::SafeDownCast(x) );
    } else {
        mapper->SetInputData( vtkUnstructuredGrid::SafeDownCast(inputGeomerty) );
    }

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0,0,0); // Background color white
    // renderer->SetBackgroundAlpha(1);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer->AddActor(actor);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize( this->Resolution );
    renderWindow->SetAlphaBitPlanes(0); //enable usage of alpha channel
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetParallelProjection(true);
    camera->SetParallelScale( this->CamScale );
    camera->SetPosition( this->CamPosition );
    camera->SetFocalPoint( this->CamFocus );
    camera->SetClippingRange( this->CamNearFar );
    renderer->SetActiveCamera(camera);
    renderWindow->Render();

    {
        stringstream msg;
        t0 = t.getElapsedTime();
        msg<<"[ttkCinemaImageExport] VTK Rendering Pipeline initialized in " << t0 << " s."<<endl;

        dMsg(cout, msg.str(), timeMsg);
    }

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput( renderWindow );
    windowToImageFilter->SetInputBufferTypeToZBuffer(); //also record the alpha (transparency) channel

    size_t n = inputGrid->GetNumberOfPoints();
    {
        auto ch = vtkSmartPointer<vtkDoubleArray>::New();
        ch->SetName("CamHeight");
        ch->SetNumberOfValues(1);
        ch->SetValue(0, this->CamScale*2.);

        auto cnf = vtkSmartPointer<vtkDoubleArray>::New();
        cnf->SetName("CamNearFar");
        cnf->SetNumberOfValues(2);
        cnf->SetValue(0, this->CamNearFar[0]);
        cnf->SetValue(1, this->CamNearFar[1]);

        auto cr = vtkSmartPointer<vtkDoubleArray>::New();
        cr->SetName("CamRes");
        cr->SetNumberOfValues(2);
        cr->SetValue(0, this->Resolution[0]);
        cr->SetValue(1, this->Resolution[1]);


        for(size_t i=0; i<n; i++){

            inputGrid->GetPoint(i, this->CamPosition);
            camera->SetPosition( this->CamPosition );

            renderWindow->Render();
            windowToImageFilter->Modified();
            windowToImageFilter->Update();


            auto outputImage = vtkSmartPointer<vtkImageData>::New();
            outputImage->DeepCopy(windowToImageFilter->GetOutput());

            auto depthValues = outputImage->GetPointData()->GetArray("ImageScalars");
            depthValues->SetName("DepthValues");

            auto outputImageFD = outputImage->GetFieldData();

            // Statics
            outputImageFD->AddArray( ch );
            outputImageFD->AddArray( cnf );
            outputImageFD->AddArray( cr );

            // Position
            auto cp = vtkSmartPointer<vtkDoubleArray>::New();
            cp->SetName("CamPosition");
            cp->SetNumberOfValues(3);
            cp->SetValue(0, this->CamPosition[0]);
            cp->SetValue(1, this->CamPosition[1]);
            cp->SetValue(2, this->CamPosition[2]);
            outputImageFD->AddArray( cp );

            // Dir
            auto cd = vtkSmartPointer<vtkDoubleArray>::New();
            cd->SetName("CamDirection");
            cd->SetNumberOfValues(3);
            double tempCD[3] = {
                this->CamFocus[0]-this->CamPosition[0],
                this->CamFocus[1]-this->CamPosition[1],
                this->CamFocus[2]-this->CamPosition[2]
            };
            vtkMath::Normalize(tempCD);
            cd->SetValue(0, tempCD[0]);
            cd->SetValue(1, tempCD[1]);
            cd->SetValue(2, tempCD[2]);
            outputImageFD->AddArray( cd );

            // Up
            auto upd = camera->GetViewUp();
            auto cu = vtkSmartPointer<vtkDoubleArray>::New();
            cu->SetName("CamUp");
            cu->SetNumberOfValues(3);
            cu->SetValue(0, upd[0]);
            cu->SetValue(1, upd[1]);
            cu->SetValue(2, upd[2]);
            outputImageFD->AddArray( cu );

            output->SetBlock(i, outputImage);
        }
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaImageExport] ------------------------------------" << endl;
        msg << "[ttkCinemaImageExport] " << n << " Images rendered" << endl;
        msg << "[ttkCinemaImageExport]   time: " << (t.getElapsedTime()-t0) << " s" << endl;
        msg << "[ttkCinemaImageExport] memory: " << m.getElapsedUsage() << " MB" << endl;
        msg << "[ttkCinemaImageExport] ------------------------------------" << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}