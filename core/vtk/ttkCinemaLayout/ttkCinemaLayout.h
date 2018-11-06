/// \ingroup vtk
/// \class ttkCinemaLayout
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.10.2018
///
/// \brief TTK VTK-filter that computes an image grid layout.
///
/// This filter computes a grid layout for images stored as blocks of a vtkMultiBlockDataSet.
///
/// \param Input vtkMultiBlockDataSet
/// \param Output vtkMultiBlockDataSet

#pragma once

// VTK includes
#include <vtkXMLPMultiBlockDataWriter.h>
#include <vtkInformation.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaLayout
#else
class ttkCinemaLayout
#endif
: public vtkXMLPMultiBlockDataWriter, public ttk::Wrapper{

    public:

        static ttkCinemaLayout* New();
        vtkTypeMacro(ttkCinemaLayout, vtkXMLPMultiBlockDataWriter)

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
            return 1;
        }

    protected:

        ttkCinemaLayout(){
            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkCinemaLayout(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

    private:

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};