/// \ingroup vtk
/// \class ttkHelloWorld
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the helloWorld processing package.
///
/// VTK wrapping code for the @HelloWorld package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::HelloWorld
#pragma once

#include <vtkInformation.h>
#include <ttkTriangulationAlgorithm.h>

// ttk includes
#include <HelloWorld.h>
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

      private:
        ttk::helloWorld::HelloWorld helloWorld_;
};
