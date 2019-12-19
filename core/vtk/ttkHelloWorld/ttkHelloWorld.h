/// \ingroup vtk
/// \class ttkHelloWorld
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2019.
///
/// \brief TTK VTK-filter that wraps the ttk::HelloWorld module.
///
/// This VTK filter uses the ttk::HelloWorld module to compute the bounding box of a vtkDataSet, which is returned as a vtkUnstructuredGrid.
///
/// \param Input vtkDataSet whose bounding box will be computed.
/// \param Output vtkUnstructuredGrid that corresponds to bounding box of the input.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::HelloWorld
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkHelloWorldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <HelloWorld.h>

class TTKHELLOWORLD_EXPORT ttkHelloWorld
    : public ttkAlgorithm    // we inherit from the generic ttkAlgorithm class
    , public ttk::HelloWorld // and we inherit from the base class
{
    private:
        /**
         * TODO 1: Add all filter parameters only as private member variables.
         */
        float Scale{1};

    public:
        /**
         * TODO 2: Automatically generate getters and setters of filter
         * parameters via vtkMacros.
         */
        vtkSetMacro(Scale, float);
        vtkGetMacro(Scale, float);

        /**
         * This static method and the macro below are VTK conventions. You don't
         * have to modify this.
         */
        static ttkHelloWorld *New();
        vtkTypeMacro(ttkHelloWorld, ttkAlgorithm);

    protected:
        /**
         * TODO 3: Implement the filter constructor and destructor (see cpp file)
         */
        ttkHelloWorld();
        ~ttkHelloWorld() override;

        /**
        * TODO 4: Specify the required input data type of each input port (see cpp file)
        */
        int FillInputPortInformation(int port, vtkInformation* info) override;

        /**
        * TODO 5: Specify the data object type of each output port (see cpp file)
        */
        int FillOutputPortInformation(int port, vtkInformation* info) override;

        /**
        * TODO 6: Implement the purpose of this filter (see cpp file)
        */
        int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector
        ) override;
};