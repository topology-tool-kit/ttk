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

vtkSmartPointer<vtkCamera> createCamera(double* camNearFar, double* camFocus, double camHeight){
    auto camera = vtkSmartPointer<vtkCamera>::New();

    camera->SetParallelProjection(true);
    camera->SetClippingRange( camNearFar );
    camera->SetFocalPoint( camFocus );
    camera->SetParallelScale( camHeight*0.5 ); // *0.5 to convert CamHeight to weird VTK convention

    return camera;
}

vtkSmartPointer<vtkRenderer> createRenderer(vtkSmartPointer<vtkActor> actor, vtkSmartPointer<vtkCamera> camera){
    auto renderer = vtkSmartPointer<vtkRenderer>::New();

    renderer->SetBackground(0,0,0); // Background color black
    renderer->AddActor(actor);
    renderer->SetActiveCamera(camera);

    return renderer;
}

vtkSmartPointer<vtkRenderWindow> createWindow(int* resolution, vtkSmartPointer<vtkRenderer> renderer){
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

    renderWindow->SetSize( resolution );
    renderWindow->SetMultiSamples( 0 ); // Disable AA
    renderWindow->AddRenderer( renderer );

    return renderWindow;
}

vector< vtkValuePass* > addValuePasses(vtkSmartPointer<vtkRenderer> renderer, vtkPolyData* poly){
    auto pointData = poly->GetPointData();
    auto cellData = poly->GetCellData();


    vector< vtkValuePass* > passesArray;

    auto passes = vtkSmartPointer<vtkRenderPassCollection>::New();

    double minmax[2];

    // Point Data
    {
        size_t n = pointData->GetNumberOfArrays();
        for(size_t i=0; i<n; i++){
            auto values = pointData->GetArray(i);

            auto valuePass = vtkSmartPointer<vtkValuePass>::New();
            valuePass->SetInputArrayToProcess(VTK_SCALAR_MODE_USE_POINT_FIELD_DATA, i);
            valuePass->SetRenderingMode(2);
            valuePass->SetInputComponentToProcess(0); // TODO Iterate over components

            values->GetRange(minmax);
            valuePass->SetScalarRange(minmax[0], minmax[1]);

            passes->AddItem(valuePass);
            passesArray.push_back( valuePass );
        }

    }

    // TODO Cell Data
    {
        size_t m = cellData->GetNumberOfArrays();
    }

    auto sequence = vtkSmartPointer<vtkSequencePass>::New();
    sequence->SetPasses(passes);

    auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
    cameraPass->SetDelegatePass(sequence);

    auto glRenderer = vtkOpenGLRenderer::SafeDownCast( renderer );
    glRenderer->SetPass(cameraPass);

    return passesArray;
}

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
    // Initialize Renderer
    // -------------------------------------------------------------------------

    // Insert InputDataObject into MultiBlockDataSet
    auto inputMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    inputMultiBlock->SetBlock(0, inputObject);

    // Create Mapper and Actor
    auto toPoly = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
    toPoly->SetInputData( inputMultiBlock );
    toPoly->Update();

    auto poly = toPoly->GetOutput();

    auto mapper = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
    mapper->SetInputConnection(toPoly->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    // Setup Camera
    auto camera = createCamera(this->CamNearFar, this->CamFocus, this->CamHeight);

    // Create Renderers and Windows
    auto depthRenderer = createRenderer( actor, camera );
    auto depthWindow = createWindow( this->Resolution, depthRenderer );

    auto valueRenderer = createRenderer( actor, camera );
    auto valueWindow = createWindow( this->Resolution, valueRenderer );
    auto valuePasses = addValuePasses(valueRenderer, poly);
    bool renderValues = valuePasses.size()>0;

    // First Render Call to Update Everything
    depthWindow->Render();

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
    windowToImageFilter->SetInput( depthWindow );
    windowToImageFilter->SetInputBufferTypeToZBuffer(); // Set output to depth buffer

    // Prepare Field Data
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
        inputGrid->GetPoint(i, camPosition);
        camera->SetPosition( camPosition );

        // TODO: In the future it shoud be possible to override focus, res,...

        // Render Image
        depthWindow->Render();

            // vtkSmartPointer<vtkWindowToImageFilter> grabber =
            // vtkSmartPointer<vtkWindowToImageFilter>::New();
            //   grabber->SetInput(renderWindow);
            //   grabber->Update();
            //   vtkImageData *id = grabber->GetOutput();
            //   //id->PrintSelf(cerr, vtkIndent(0));

            //   vtkUnsignedCharArray *ar =
            //     vtkArrayDownCast<vtkUnsignedCharArray>(id->GetPointData()->GetArray("ImageScalars"));
            //   unsigned char *ptr = static_cast<unsigned char*>(ar->GetVoidPointer(0));
            //   double value;
            //   for (int i = 0; i < id->GetNumberOfPoints(); i++)
            //   {
            //     valuePass->ColorToValue(ptr, 0, 1, value);
            //         ptr+=3;
            //   }

        windowToImageFilter->Modified();
        windowToImageFilter->Update();

        // Get Image
        auto outputImage = vtkSmartPointer<vtkImageData>::New();
        outputImage->DeepCopy(windowToImageFilter->GetOutput());

        // Rename Scalar Values
        auto depthValues = outputImage->GetPointData()->GetArray("ImageScalars");
        depthValues->SetName("DepthValues");

        // Additional Point Data
            valueWindow->Render();
        // if(renderValues){
        //     // for(auto pass: valuePasses){
        //     // }
        // }

        // auto outputImagePD = outputImage->GetPointData();
        // vtkSmartPointer<vtkFloatArray> ele = vtkSmartPointer<vtkFloatArray>::New();
        // ele->DeepCopy( valuePass->GetFloatImageDataArray( renderer ) );
        // ele->SetName("xxx");
        // outputImagePD->AddArray( ele );

        // Set Field Data of Output Image
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
        auto upd = camera->GetViewUp();
        auto cu = vtkSmartPointer<vtkDoubleArray>::New();
        cu->SetName("CamUp");
        cu->SetNumberOfValues(3);
        cu->SetValue(0, upd[0]);
        cu->SetValue(1, upd[1]);
        cu->SetValue(2, upd[2]);
        outputImageFD->AddArray( cu );

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