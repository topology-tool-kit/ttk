#include <ttkCinemaImaging.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkActor.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkCamera.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <vtkFloatArray.h>
#include <vtkValuePass.h>
#include <vtkRenderPassCollection.h>
#include <vtkCameraPass.h>
#include <vtkSequencePass.h>
#include <vtkOpenGLRenderer.h>

#include <set>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaImaging)

int ttkCinemaImaging::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaImaging] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;
    Timer t;
    double t0=0;

    // Get Input / Output
    vtkInformation* inputGeomertyInfo = inputVector[0]->GetInformationObject(0);
    auto inputObject = inputGeomertyInfo->Get(vtkDataObject::DATA_OBJECT());

    vtkInformation* inGridInfo = inputVector[1]->GetInformationObject(0);
    auto inputGrid = vtkPointSet::SafeDownCast( inGridInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // -------------------------------------------------------------------------
    // Initialize Render Object
    // -------------------------------------------------------------------------

    // Insert InputDataObject into MultiBlockDataSet
    auto inputMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    inputMultiBlock->SetBlock(0, inputObject);

    // Convert MultiBlock to PolyData
    auto toPoly = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
    toPoly->SetInputData( inputMultiBlock );
    toPoly->Update();
    auto poly = toPoly->GetOutput();

    // -------------------------------------------------------------------------
    // Initialize Depth Renderer
    // -------------------------------------------------------------------------

    // Mapper
    auto mapper0 = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
    mapper0->SetInputConnection(toPoly->GetOutputPort());

    // Actor
    auto actor0 = vtkSmartPointer<vtkActor>::New();
    actor0->SetMapper(mapper0);

    // Camera
    auto camera0 = vtkSmartPointer<vtkCamera>::New();
    camera0->SetParallelProjection(true);
    camera0->SetClippingRange( this->CamNearFar );
    camera0->SetFocalPoint( this->CamFocus );
    camera0->SetParallelScale( this->CamHeight*0.5 ); // *0.5 to convert CamHeight to weird VTK convention

    // Renderer
    auto renderer0 = vtkSmartPointer<vtkRenderer>::New();
    renderer0->SetBackground(0,0,0);
    renderer0->AddActor(actor0);
    renderer0->SetActiveCamera(camera0);

    // Window
    auto renderWindow0 = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow0->SetSize( this->Resolution );
    renderWindow0->SetMultiSamples( 0 ); // Disable AA
    renderWindow0->AddRenderer( renderer0 );

    // -------------------------------------------------------------------------
    // Initialize Value Renderer
    // -------------------------------------------------------------------------

    // Mapper
    auto mapper1 = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
    mapper1->SetInputConnection(toPoly->GetOutputPort());

    // Actor
    auto actor1 = vtkSmartPointer<vtkActor>::New();
    actor1->SetMapper(mapper1);

    // Camera
    auto camera1 = vtkSmartPointer<vtkCamera>::New();
    camera1->SetParallelProjection(true);
    camera1->SetClippingRange( this->CamNearFar );
    camera1->SetFocalPoint( this->CamFocus );
    camera1->SetParallelScale( this->CamHeight*0.5 ); // *0.5 to convert CamHeight to weird VTK convention

    // Renderer
    auto renderer1 = vtkSmartPointer<vtkRenderer>::New();
    renderer1->SetBackground(0,0,0);
    renderer1->AddActor(actor1);
    renderer1->SetActiveCamera(camera1);

    // Window
    auto renderWindow1 = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow1->SetSize( this->Resolution );
    renderWindow1->SetMultiSamples( 0 ); // Disable AA
    renderWindow1->AddRenderer( renderer1 );

    // Value Passes
    vector< pair<vtkValuePass*,string> > valuePasses;
    {
        auto passes = vtkSmartPointer<vtkRenderPassCollection>::New();

        auto pointData = poly->GetPointData();
        size_t n = pointData->GetNumberOfArrays();

        double minmax[2];
        for(size_t i=0; i<n; i++){
            auto values = pointData->GetArray(i);
            values->GetRange(minmax);

            auto valuePass = vtkSmartPointer<vtkValuePass>::New();
            valuePass->SetInputArrayToProcess(VTK_SCALAR_MODE_USE_POINT_FIELD_DATA, i);
            valuePass->SetRenderingMode(2);
            valuePass->SetInputComponentToProcess(0); // TODO Iterate components
            valuePass->SetScalarRange(minmax[0], minmax[1]);

            passes->AddItem(valuePass);
            valuePasses.push_back( make_pair(valuePass, values->GetName()) );
        }

        // TODO Iterate Cell Data and Components

        auto sequence = vtkSmartPointer<vtkSequencePass>::New();
        sequence->SetPasses(passes);

        auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
        cameraPass->SetDelegatePass(sequence);

        auto glRenderer = vtkOpenGLRenderer::SafeDownCast( renderer1 );
        glRenderer->SetPass(cameraPass);
    }
    bool renderValuePasses = valuePasses.size()>0;

    renderWindow0->Render();
    if(renderValuePasses) renderWindow1->Render();

    // Print Status
    {
        stringstream msg;
        t0 = t.getElapsedTime();
        msg<<"[ttkCinemaImaging] VTK Rendering Pipeline initialized in " << t0 << " s."<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // -------------------------------------------------------------------------
    // Render Images for all Camera Locations
    // -------------------------------------------------------------------------

    // Create vtkWindowToImageFilter
    auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput( renderWindow0 );
    windowToImageFilter->SetInputBufferTypeToZBuffer(); // Set output to depth buffer

    // Prepare Field Data for Depth Values
    auto ch = vtkSmartPointer<vtkDoubleArray>::New();
    ch->SetName("CamHeight");
    ch->SetNumberOfValues(1);
    ch->SetValue(0, this->CamHeight);

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

    // Iterate over Locations
    double camPosition[3] = {0,0,0};
    size_t n = inputGrid->GetNumberOfPoints();

    for(size_t i=0; i<n; i++){
        // Set Camera Position
        // TODO: In the future it shoud be possible to override focus, res,...
        inputGrid->GetPoint(i, camPosition);
        camera0->SetPosition( camPosition );
        camera1->SetPosition( camPosition );

        // Render Depth Image
        renderWindow0->Render();

        // Initialize Output Image
        auto outputImage = vtkSmartPointer<vtkImageData>::New();
        {
            windowToImageFilter->Modified();
            windowToImageFilter->Update();
            outputImage->DeepCopy(windowToImageFilter->GetOutput());
            auto depthValues = outputImage->GetPointData()->GetAbstractArray("ImageScalars");
            depthValues->SetName("Depth");
        }

        // Add Field Data
        {
            auto outputImageFD = outputImage->GetFieldData();

            // Camera Parameters
            outputImageFD->AddArray( ch );
            outputImageFD->AddArray( cnf );
            outputImageFD->AddArray( cr );

            // Position
            auto cp = vtkSmartPointer<vtkDoubleArray>::New();
            cp->SetName("CamPosition");
            cp->SetNumberOfValues(3);
            cp->SetValue(0, camPosition[0]);
            cp->SetValue(1, camPosition[1]);
            cp->SetValue(2, camPosition[2]);
            outputImageFD->AddArray( cp );

            // Dir
            auto cd = vtkSmartPointer<vtkDoubleArray>::New();
            cd->SetName("CamDirection");
            cd->SetNumberOfValues(3);
            double tempCD[3] = {
                this->CamFocus[0]-camPosition[0],
                this->CamFocus[1]-camPosition[1],
                this->CamFocus[2]-camPosition[2]
            };
            vtkMath::Normalize(tempCD);
            cd->SetValue(0, tempCD[0]);
            cd->SetValue(1, tempCD[1]);
            cd->SetValue(2, tempCD[2]);
            outputImageFD->AddArray( cd );

            // Up
            auto upd = camera0->GetViewUp();
            auto cu = vtkSmartPointer<vtkDoubleArray>::New();
            cu->SetName("CamUp");
            cu->SetNumberOfValues(3);
            cu->SetValue(0, upd[0]);
            cu->SetValue(1, upd[1]);
            cu->SetValue(2, upd[2]);
            outputImageFD->AddArray( cu );
        }

        // Add Point Data
        if(renderValuePasses){

            // Render Value Passes
            renderWindow1->Render();

            auto outputImagePD = outputImage->GetPointData();
            for(auto& passData: valuePasses){
                vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
                data->DeepCopy( passData.first->GetFloatImageDataArray( renderer1 ) );
                data->SetName( passData.second.data() );
                outputImagePD->AddArray( data );
            }
        }

        // Add Image to MultiBlock
        output->SetBlock(i, outputImage);
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaImaging] ------------------------------------" << endl;
        msg << "[ttkCinemaImaging] " << n << " Images rendered" << endl;
        msg << "[ttkCinemaImaging]   time: " << (t.getElapsedTime()-t0) << " s" << endl;
        msg << "[ttkCinemaImaging] memory: " << m.getElapsedUsage() << " MB" << endl;
        msg << "[ttkCinemaImaging] ------------------------------------" << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}