/// \ingroup vtk
/// \class ttkEndFor
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that requests more for executions
///
/// This filter requests more for executions as long as the maximum number of elements is not reached. This filter works in conjunction with  the ttkForEachRow or ttkForEachBlock filter.
///
/// \param Input vtkDataSet Data that will be pushed trough after all iterations.
/// \param Output vtkMultiBlockDataSet Shallow copy of the input

#pragma once

// VTK includes
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkEndFor
#else
class ttkEndFor
#endif
: public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

    public:

        static ttkEndFor* New();
        vtkTypeMacro(ttkEndFor, vtkMultiBlockDataSetAlgorithm)

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

        int FillInputPortInformation(int port, vtkInformation* info) override {
            switch(port)
                case 1: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation* info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
            return 1;
        }

    protected:

        ttkEndFor(){
            nextIndex = 0;

            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkEndFor(){};

        bool UseAllCores;
        int ThreadNumber;

        double nextIndex;

        int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;
        int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;
        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};