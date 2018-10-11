/// \ingroup vtk
/// \class ttkDepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.07.2018
///
/// \brief TTK VTK-filter that approximates the geomerty that is depicted by a set of depth images.
///
/// VTK wrapping code for the @DepthImageBasedGeometryApproximation package.
///
/// \param Input Depth Images (vtkMultiBlockDataSet)
/// \param Output A set of unstructured grids where each grid corresponds to a depth image (vtkMultiBlockDataSet)
///
/// \sa ttk::DepthImageBasedGeometryApproximation
#pragma once

#include                  <vtkMultiBlockDataSet.h>
#include                  <vtkMultiBlockDataSetAlgorithm.h>
#include                  <vtkInformation.h>

#include                  <DepthImageBasedGeometryApproximation.h>
#include                  <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDepthImageBasedGeometryApproximation
#else
class ttkDepthImageBasedGeometryApproximation
#endif
: public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

    public:

        static ttkDepthImageBasedGeometryApproximation* New();
        vtkTypeMacro(ttkDepthImageBasedGeometryApproximation, vtkMultiBlockDataSetAlgorithm)

        vtkSetMacro(SubSampling, int);
        vtkGetMacro(SubSampling, int);

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

        ttkDepthImageBasedGeometryApproximation(){
            SubSampling = 0;

            UseAllCores = false;
            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkDepthImageBasedGeometryApproximation(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    private:

        int SubSampling;
        ttk::DepthImageBasedGeometryApproximation depthImageBasedGeometryApproximation_;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
