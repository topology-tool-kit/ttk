#pragma once

#include <vtkInformation.h>
#include <ttkTriangulationFilter.h>

// VTK Module
#include <ttkHelloWorldModule.h>

class TTKHELLOWORLD_EXPORT ttkHelloWorld
    : public ttkTriangulationFilter
{

    public:
        static ttkHelloWorld *New();
        vtkTypeMacro(ttkHelloWorld, ttkTriangulationFilter);

    protected:
        ttkHelloWorld();
        ~ttkHelloWorld();

        int FillInputPortInformation(int port, vtkInformation* info) override;
        int FillOutputPortInformation(int port, vtkInformation* info) override;

        int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector
        ) override;
};
