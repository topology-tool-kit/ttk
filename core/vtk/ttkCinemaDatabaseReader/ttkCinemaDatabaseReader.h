/// \ingroup vtk
/// \class ttkCinemaDatabaseReader
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database
///
/// \param Output content of the data.csv file of the database in form of a vtkTable
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \sa ttk::CinemaDatabaseReader
#pragma once

// VTK includes
#include <vtkTableReader.h>
#include <vtkInformation.h>
#include <vtkTable.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaDatabaseReader
#else
class ttkCinemaDatabaseReader
#endif
: public vtkTableReader, public ttk::Wrapper{

    public:

        static ttkCinemaDatabaseReader* New();
        vtkTypeMacro(ttkCinemaDatabaseReader, vtkTableReader)

        vtkSetMacro(DatabasePath, std::string);
        vtkGetMacro(DatabasePath, std::string);

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

        ttkCinemaDatabaseReader(){
            DatabasePath = "";

            UseAllCores = true;
            ThreadNumber = 1;

            SetNumberOfInputPorts(0);
            SetNumberOfOutputPorts(1);
        }
        ~ttkCinemaDatabaseReader(){};

        bool UseAllCores;
        int ThreadNumber;

        int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    private:

        std::string DatabasePath;

        bool needsToAbort() override { return GetAbortExecute();};
        int updateProgress(const float &progress) override {
            UpdateProgress(progress);
            return 0;
        };
};
