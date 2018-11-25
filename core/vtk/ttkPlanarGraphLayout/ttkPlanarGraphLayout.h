/// \ingroup vtk
/// \class ttkPlanarGraphLayout
/// \author Wiebke Koepp (wiebke.koepp@gmail.com) and Jonas Lukasczyk (jl@jluk.de)
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that computes a planar graph layout.
///
/// VTK wrapping code for the @PlanarGraphLayout package.
///
/// This filter computes a planar graph layout of a vtkUnstructuredGrid.
///
/// \param Input Graph. (vtkUnstructuredGrid)
/// \param Output Graph (vtkUnstructuredGrid)
///
/// \sa ttk::PlanarGraphLayout

#pragma once

// VTK includes
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <PlanarGraphLayout.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPlanarGraphLayout
#else
class ttkPlanarGraphLayout
#endif
: public vtkUnstructuredGridAlgorithm, public ttk::Wrapper{

    public:

        static ttkPlanarGraphLayout* New();
        vtkTypeMacro(ttkPlanarGraphLayout, vtkUnstructuredGridAlgorithm)

        vtkSetMacro(AxisFieldName, std::string);
        vtkGetMacro(AxisFieldName, std::string);

        vtkSetMacro(UseSizes, bool);
        vtkGetMacro(UseSizes, bool);
        vtkSetMacro(SizeFieldName, std::string);
        vtkGetMacro(SizeFieldName, std::string);

        vtkSetMacro(UseLevels, bool);
        vtkGetMacro(UseLevels, bool);
        vtkSetMacro(LevelFieldName, std::string);
        vtkGetMacro(LevelFieldName, std::string);

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
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port)
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
            return 1;
        }

    protected:

        ttkPlanarGraphLayout(){
            SetAxisFieldName("TimeIndex");
            SetUseSizes(false);
            SetSizeFieldName("Size");
            SetUseLevels(false);
            SetLevelFieldName("LevelIndex");

            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkPlanarGraphLayout(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        std::string AxisFieldName;

        bool UseSizes;
        std::string SizeFieldName;

        bool UseLevels;
        std::string LevelFieldName;

        ttk::PlanarGraphLayout planarGraphLayout_;

        bool needsToAbort() override { return GetAbortExecute(); };
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
