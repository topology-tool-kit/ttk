/// \ingroup vtk
/// \class ttkPlanarGraphLayout
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK VTK-filter that computes a planar graph layout.
///
/// This filter computes a planar graph layout of a 'vtkUnstructuredGrid'. To improve the quality of the layout it is possible to pass additional field data to the algorithm:
///
/// 1) Sequences: Points are positioned along the x-axis based on a sequence (e.g., time indicies or scalar values).
/// 2) Sizes: Points cover space on the y-axis based on their size.
/// 3) Branches: Points with the same branch label are positioned on straight lines.
/// 4) Levels: The layout of points with the same level label are computed individually and afterwards nested based on the level hierarchy. This makes it possible to draw nested graphs where each level is a layer of the resulting graph.
///
/// VTK wrapping code for the @PlanarGraphLayout package.
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

        // getters and setters for optional field data
        vtkSetMacro(UseSequences, bool);
        vtkGetMacro(UseSequences, bool);
        vtkSetMacro(SequenceFieldName, std::string);
        vtkGetMacro(SequenceFieldName, std::string);

        vtkSetMacro(UseSizes, bool);
        vtkGetMacro(UseSizes, bool);
        vtkSetMacro(SizeFieldName, std::string);
        vtkGetMacro(SizeFieldName, std::string);

        vtkSetMacro(UseBranches, bool);
        vtkGetMacro(UseBranches, bool);
        vtkSetMacro(BranchFieldName, std::string);
        vtkGetMacro(BranchFieldName, std::string);

        vtkSetMacro(UseLevels, bool);
        vtkGetMacro(UseLevels, bool);
        vtkSetMacro(LevelFieldName, std::string);
        vtkGetMacro(LevelFieldName, std::string);

        // getters and setters for output field name
        vtkSetMacro(OutputFieldName, std::string);
        vtkGetMacro(OutputFieldName, std::string);

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
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); break;
                default: return 0;
            }
            return 1;
        }

        int FillOutputPortInformation(int port, vtkInformation *info) override {
            switch(port){
                case 0: info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); break;
                default: return 0;
            }
            return 1;
        }

    protected:

        ttkPlanarGraphLayout(){
            SetUseSequences(false);
            SetSequenceFieldName("");
            SetUseSizes(false);
            SetSizeFieldName("");
            SetUseBranches(false);
            SetBranchFieldName("");
            SetUseLevels(false);
            SetLevelFieldName("");

            SetOutputFieldName("Layout");

            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkPlanarGraphLayout(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        // optional field data
        bool UseSequences;
        std::string SequenceFieldName;
        bool UseSizes;
        std::string SizeFieldName;
        bool UseBranches;
        std::string BranchFieldName;
        bool UseLevels;
        std::string LevelFieldName;

        // output field name
        std::string OutputFieldName;

        // base code
        ttk::PlanarGraphLayout planarGraphLayout;

        bool needsToAbort() override { return GetAbortExecute(); };
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
