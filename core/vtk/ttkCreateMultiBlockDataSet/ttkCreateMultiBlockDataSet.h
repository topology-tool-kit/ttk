/// \ingroup vtk
/// \class ttkCreateMultiBlockDataSet
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.9.2018
///
/// \brief TTK VTK-filter that generates images of a vtkDataSet.
///
/// This filter takes images of a vtkDataObject from positions specified on a vtkPointSet. Each image will be a block of a vtkMultiBlockDataSet where block order corresponds to point order. Each sample point can optionally have vtkDoubleArrays to override the default rendering parameters, i.e, the resolution, focus, clipping planes, and viewport height.
///
/// VTK wrapping code for the @CreateMultiBlockDataSet package.
///
/// \param Input vtkDataObject that will be depicted (vtkDataObject)
/// \param Input vtkPointSet that records the camera sampling locations (vtkPointSet)
/// \param Output vtkMultiBlockDataSet that represents a list of images (vtkMultiBlockDataSet)

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCreateMultiBlockDataSet
#else
class ttkCreateMultiBlockDataSet
#endif
: public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

    public:

        static ttkCreateMultiBlockDataSet* New();
        vtkTypeMacro(ttkCreateMultiBlockDataSet, vtkMultiBlockDataSetAlgorithm)

        // default ttk setters
        vtkSetMacro(debugLevel_, int);
        void SetThreads(){
            threadNumber_ = !UseAllCores ? ThreadNumber : ttk::OsCall::getNumberOfCores();
            Modified();
        }
        void SetThreadNumber(int threadNumber){
            ThreadNumber = threadNumber;
            SetThreads();
        }
        void SetUseAllCores(bool onOff){
            UseAllCores = onOff;
            SetThreads();
        }
        // end of default ttk setters

        int FillInputPortInformation(int port, vtkInformation *info) override {
            // All ports have the same generic type and are optional
            info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
            info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1 );
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
            return 1;
        }

    protected:

        ttkCreateMultiBlockDataSet(){
            UseAllCores = false;

            SetNumberOfInputPorts(5);
            SetNumberOfOutputPorts(1);
        }
        ~ttkCreateMultiBlockDataSet(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    private:

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
