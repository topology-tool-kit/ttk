/// \ingroup vtk
/// \class ttkCinemaQueryReader
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that wraps the cinemaQueryReader processing package.
///
/// VTK wrapping code for the @CinemaQueryReader package.
///
/// param Input A table where a column contains filepaths (vtkTable)
/// param Output A vtkMultiBlockDataSet where each block represents a row of the input table (vtkMultiBlockDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// sa ttk::CinemaQueryReader
#pragma once

#include                  <vtkTable.h>
#include                  <vtkImageData.h>

#include                  <vtkMultiBlockDataSetAlgorithm.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkInformation.h>
#include                  <vtkMultiBlockDataSet.h>

// ttk code includes
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkImageAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaQueryReader
#else
class ttkCinemaQueryReader
#endif
  : public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkCinemaQueryReader* New();
    vtkTypeMacro(ttkCinemaQueryReader, vtkMultiBlockDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else
        threadNumber_ = ttk::OsCall::getNumberOfCores();
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

    // vtkGetMacro(FilepathColumnName, std::string);
    // vtkSetMacro(FilepathColumnName, std::string);

    void SetFilepathColumnName(int idx, int port, int connection, int fieldAssociation, const char *name){
        this->FilepathColumnName = std::string(name);
        this->Modified();
    };


    int FillInputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
          break;
        default:
          break;
      }
      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
          break;
        default:
          break;
      }
      return 1;
    }

  protected:

    ttkCinemaQueryReader(){

        UseAllCores = false;

        SetNumberOfInputPorts(1);
        SetNumberOfOutputPorts(1);
    }
    ~ttkCinemaQueryReader(){};

    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    bool UseAllCores;
    int ThreadNumber;

    std::string FilepathColumnName;

  private:

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };

};
