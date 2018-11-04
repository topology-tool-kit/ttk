/// \ingroup vtk
/// \class ttkAddFieldData
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads the data products that are referenced in a vtkTable.
///
/// This filter reads the products that are referenced in a vtkTable. The results are stored in a vtkMultiBlockDataSet where each block corresponds to a row of the table with consistent ordering.
///
/// \param Input vtkTable that contains data product references (vtkTable)
/// \param Output vtkMultiBlockDataSet where each block is a referenced product of an input table row (vtkMultiBlockDataSet)

#pragma once

// VTK includes
#include <vtkTableAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkAddFieldData
#else
class ttkAddFieldData
#endif
: public vtkTableAlgorithm, public ttk::Wrapper{

    public:

        static ttkAddFieldData* New();
        vtkTypeMacro(ttkAddFieldData, vtkTableAlgorithm)

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
            switch(port){
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject"); break;
                case 1: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject"); break;
            }
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
            return 1;
        }

    protected:

        ttkAddFieldData(){
            UseAllCores = false;

            SetNumberOfInputPorts(2);
            SetNumberOfOutputPorts(1);
        }
        ~ttkAddFieldData(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        std::string FieldDataString;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};