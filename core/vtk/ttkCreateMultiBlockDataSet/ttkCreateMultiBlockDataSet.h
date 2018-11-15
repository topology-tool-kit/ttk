/// \ingroup vtk
/// \class ttkCreateMultiBlockDataSet
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that stores up to 5 vtkDataObjects as blocks of a vtkMultiBlockDataSet.
///
/// This filter takes up to 5 vtkDataObjects and stores them as blocks of a newly created vtkMultiBlockDataSet.
///
/// \param Input vtkDataObject
/// \param Input vtkDataObject
/// \param Input vtkDataObject
/// \param Input vtkDataObject
/// \param Input vtkDataObject
/// \param Output vtkMultiBlockDataSet

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
            if(port>4)
                return 0;
            // All ports have the same generic type and are optional
            info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
            info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1 );
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port){
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet"); break;
                default: return 0;
            }
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

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        bool needsToAbort() override { return GetAbortExecute(); };
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
