/// \ingroup vtk
/// \class ttkMeshGraph
/// \author Wiebke Koepp (wiebke.koepp@gmail.com) and Jonas Lukasczyk (jl@jluk.de)
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that TODO.
///
/// VTK wrapping code for the @MeshGraph package.
///
/// This filter TODO
///
/// \param Input Graph. (vtkUnstructuredGrid)
/// \param Output Graph (vtkUnstructuredGrid)
///
/// \sa ttk::MeshGraph

#pragma once

// VTK includes
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <MeshGraph.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMeshGraph
#else
class ttkMeshGraph
#endif
: public vtkUnstructuredGridAlgorithm, public ttk::Wrapper{

    public:

        static ttkMeshGraph* New();
        vtkTypeMacro(ttkMeshGraph, vtkUnstructuredGridAlgorithm)

        vtkSetMacro(SizeFieldName, std::string);
        vtkGetMacro(SizeFieldName, std::string);

        vtkSetMacro(PrimaryAxis, int);
        vtkGetMacro(PrimaryAxis, int);

        vtkSetMacro(SecondaryAxis, int);
        vtkGetMacro(SecondaryAxis, int);

        vtkSetMacro(Subdivisions, int);
        vtkGetMacro(Subdivisions, int);

        vtkSetMacro(SizeScale, float);
        vtkGetMacro(SizeScale, float);

        vtkSetMacro(Tetrahedralize, bool);
        vtkGetMacro(Tetrahedralize, bool);

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

        ttkMeshGraph(){
            SetSizeFieldName("Size");
            SetPrimaryAxis(0);
            SetSecondaryAxis(1);
            SetSubdivisions(0);
            SetSizeScale(1);

            UseAllCores = false;

            SetNumberOfInputPorts(1);
            SetNumberOfOutputPorts(1);
        }
        ~ttkMeshGraph(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

    private:

        std::string SizeFieldName;
        int PrimaryAxis;
        int SecondaryAxis;
        int Subdivisions;
        float SizeScale;
        bool Tetrahedralize;

        ttk::MeshGraph meshGraph;

        bool needsToAbort() override { return GetAbortExecute(); };
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
