/// \ingroup vtk
/// \class ttkTableSource
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.9.2018
///
/// \brief TTK VTK-filter that .
///
/// TODO
///
/// VTK wrapping code for the @TableSource package.
///
/// \param Output vtkTable

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkTableAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTableSource
#else
class ttkTableSource
#endif
: public vtkTableAlgorithm, public ttk::Wrapper{

    public:

        static ttkTableSource* New();
        vtkTypeMacro(ttkTableSource, vtkTableAlgorithm)

        vtkSetMacro(FieldDataString, std::string);
        vtkGetMacro(FieldDataString, std::string);

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
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
            return 1;
        }

    protected:

        ttkTableSource(){
            UseAllCores = false;

            SetNumberOfInputPorts(0);
            SetNumberOfOutputPorts(1);
        }
        ~ttkTableSource(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    private:

        std::string FieldDataString;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
