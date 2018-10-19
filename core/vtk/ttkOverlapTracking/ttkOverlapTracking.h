/// \ingroup vtk
/// \class ttkOverlapTracking
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that TODO
///
/// This filter TODO
///
/// VTK wrapping code for the @OverlapTracking package.
///
/// \param Input TODO
/// \param Output TODO
///
/// sa ttk::OverlapTracking

#pragma once

// VTK includes
#include <vtkTableAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <OverlapTracking.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkOverlapTracking
#else
class ttkOverlapTracking
#endif
: public vtkTableAlgorithm, public ttk::Wrapper{

    public:
        static ttkOverlapTracking* New();
        vtkTypeMacro(ttkOverlapTracking, vtkTableAlgorithm)

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
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
                // case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
            return 1;
        }

    protected:

        ttkOverlapTracking(){
            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkOverlapTracking(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    private:

        ttk::OverlapTracking overlapTracking;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
