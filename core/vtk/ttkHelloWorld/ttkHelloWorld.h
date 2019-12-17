#pragma once

#include <vtkInformation.h>
#include <ttkTriangulationAlgorithm.h>

// VTK Module
#include <ttkHelloWorldModule.h>

class TTKHELLOWORLD_EXPORT ttkHelloWorld
    : public ttkTriangulationAlgorithm
{

    public:
        static ttkHelloWorld *New();
        vtkTypeMacro(ttkHelloWorld, ttkTriangulationAlgorithm);

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
