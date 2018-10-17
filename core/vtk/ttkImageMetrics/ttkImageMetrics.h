/// \ingroup vtk
/// \class ttkImageMetrics
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.10.2018
///
/// \brief TTK VTK-filter that computes image metrics for two lists of images.
///
/// This filter takes two lists of images ("Images A" and "Images B") and computes for each pair a set of user defined image metrics.
///
/// VTK wrapping code for the @ImageMetrics package.
///
/// param Input ImagesA (vtkMultiBlockDataSet)
/// param Input ImagesB (vtkMultiBlockDataSet)
/// param Output metrics (vtkTable)
///
/// sa ttk::ImageMetrics
#pragma once

#include <vtkTableAlgorithm.h>
#include <vtkInformation.h>

#include <ImageMetrics.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkImageMetrics
#else
class ttkImageMetrics
#endif
: public vtkTableAlgorithm, public ttk::Wrapper{

    public:
        static ttkImageMetrics* New();
        vtkTypeMacro(ttkImageMetrics, vtkTableAlgorithm)

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
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");break;
                case 1: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");break;
            }
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
            return 1;
        }

    protected:

        ttkImageMetrics(){
            UseAllCores = false;

            SetNumberOfInputPorts(2);
            SetNumberOfOutputPorts(1);
        }
        ~ttkImageMetrics(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector);

    private:

        ttk::ImageMetrics imageMetrics;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
