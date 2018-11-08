/// \ingroup vtk
/// \class ttkBlockAggregator
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that iteratively appends its input to a vtkMultiBlockDataSet.
///
/// This filter iteratively appends its input as blocks to a vtkMultiBlockDataSet.
///
/// \param Input vtkDataObject that will be added as a block (vtkDataObject).
/// \param Output vtkMultiBlockDataSet containing all added blocks (vtkMultiBlockDataSet).

#pragma once

// VTK includes
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBlockAggregator
#else
class ttkBlockAggregator
#endif
: public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

    public:

        static ttkBlockAggregator* New();
        vtkTypeMacro(ttkBlockAggregator, vtkMultiBlockDataSetAlgorithm)

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
                case 1: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation* info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
            return 1;
        }

    protected:

        ttkBlockAggregator(){
            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkBlockAggregator(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        vtkSmartPointer<vtkMultiBlockDataSet> AggregatedMultiBlockDataSet;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};