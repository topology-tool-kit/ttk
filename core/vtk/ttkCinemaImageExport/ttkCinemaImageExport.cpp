#include                  <ttkCinemaImageExport.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkDataSet.h>
#include                  <vtkMultiBlockDataSet.h>
#include                  <vtkWindowToImageFilter.h>
#include                  <vtkCompositePolyDataMapper.h>
#include                  <vtkActor.h>
#include                  <vtkRenderer.h>
#include                  <vtkRenderWindow.h>
#include                  <vtkCamera.h>

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

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet* inputGeomerty = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Visualize
    vtkSmartPointer<vtkCompositePolyDataMapper> mapper = vtkSmartPointer<vtkCompositePolyDataMapper>::New();
    mapper->SetInputDataObject( inputGeomerty );

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(256,256);
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetPosition(0, 0, 20);
    camera->SetFocalPoint(0, 0, 0);
    renderer->SetActiveCamera(camera);

    // renderWindow->SetAlphaBitPlanes(0); //enable usage of alpha channel

    // vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    // vtkSmartPointer<vtkRenderWindowInteractor>::New();
    // renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->SetBackground(1,1,1); // Background color white

    renderWindow->Render();

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput( renderWindow );
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
    windowToImageFilter->Update();

    output->SetBlock(0, windowToImageFilter->GetOutput());

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaImageExport] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}