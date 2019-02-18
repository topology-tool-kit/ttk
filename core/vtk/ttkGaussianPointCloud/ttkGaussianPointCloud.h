/// \ingroup vtk
/// \class ttkGaussianPointCloud
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.10.2018
///
/// \brief TTK VTK-filter that generates an Icosphere.
///
/// VTK wrapping code for the @GaussianPointCloud package.
///
/// \sa ttk::GaussianPointCloud

#pragma once

// VTK includes
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <GaussianPointCloud.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkGaussianPointCloud
#else
class ttkGaussianPointCloud
#endif
: public vtkUnstructuredGridAlgorithm, public ttk::Wrapper{

    public:

        static ttkGaussianPointCloud* New();
        vtkTypeMacro(ttkGaussianPointCloud, vtkUnstructuredGridAlgorithm)

        vtkSetMacro(Subdivisions, int);
        vtkGetMacro(Subdivisions, int);

        vtkSetVector3Macro(Center, float);
        vtkGetVector3Macro(Center, float);

        vtkSetMacro(Radius, float);
        vtkGetMacro(Radius, float);

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
            return 0;
        }

        int FillOutputPortInformation(int port, vtkInformation* info) override {
            switch(port){
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); break;
                default: return 0;
            }
            return 1;
        }

    protected:

        ttkGaussianPointCloud(){
            SetSubdivisions( 0 );
            float center[3] = {0,0,0};
            SetCenter( center );
            SetRadius( 1 );

            UseAllCores = false;
            SetNumberOfInputPorts(0);
            SetNumberOfOutputPorts(1);
        }
        ~ttkGaussianPointCloud(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        int Subdivisions;
        float Center[3];
        float Radius;
        ttk::GaussianPointCloud gaussianPointCloud_;

        bool needsToAbort() override { return GetAbortExecute(); };
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
