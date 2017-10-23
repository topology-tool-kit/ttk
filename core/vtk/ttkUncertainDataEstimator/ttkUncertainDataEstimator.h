/// \ingroup vtkWrappers
/// \class ttkUncertainDataEstimator
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK VTK-filter that takes an input ensemble data set 
/// (represented by a list of scalar fields) and which computes various 
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \param Input0 Input ensemble scalar field #0 (vtkDataSet) 
/// \param Input1 Input ensemble scalar field #1 (vtkDataSet)\n 
/// ...\n
/// \param InputN Input ensemble scalar field #N (vtkDataSet)
/// \param Output0 Lower and upper bound fields (vtkDataSet)
/// \param Output1 Histogram estimations of the vertex probability density 
/// functions (vtkDataSet)
/// \param Output2 Mean field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example 
/// within a VTK pipeline.
///
/// \sa vtkMandatoryCriticalPoints
/// \sa ttk::UncertainDataEstimator
#ifndef _TTK_UNCERTAINDATAESTIMATOR_H
#define _TTK_UNCERTAINDATAESTIMATOR_H

// ttk code includes
#include                  <UncertainDataEstimator.h>
#include                  <Wrapper.h>

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkInformationVector.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkTable.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class VTKFILTERSCORE_EXPORT ttkUncertainDataEstimator
  : public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkUncertainDataEstimator* New();

    vtkTypeMacro(ttkUncertainDataEstimator, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
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

    // set-getters macros to define from each variable you want to access from
    // the outside (in particular from paraview) - to adapt.

    void ComputeLowerBound(bool state){
      computeLowerBound_ = state;
    }

    void ComputeUpperBound(bool state){
      computeUpperBound_ = state;
    }

    void BoundToCompute(int value){
      switch(value){
        case 0:
          ComputeLowerBound(true);
          ComputeUpperBound(true);
          break;
        case 1:
          ComputeLowerBound(true);
          ComputeUpperBound(false);
          break;
        case 2:
          ComputeLowerBound(false);
          ComputeUpperBound(true);
          break;
      }
      Modified();
    }

    void BinCount(int binCount){
      binCount_ = binCount;
      Modified();
    }

    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);


  protected:

    ttkUncertainDataEstimator();

    ~ttkUncertainDataEstimator();

    int FillInputPortInformation(int port, vtkInformation *info);
    int FillOutputPortInformation(int port, vtkInformation *info);

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);


  private:

    bool                  UseAllCores;
    int                   ThreadNumber;
    string                ScalarField;
    bool                  computeLowerBound_;
    bool                  computeUpperBound_;
    int                   binCount_;
    int                   allocatedBinCount_;
    vtkDataArray          *outputLowerBoundScalarField_;
    vtkDataArray          *outputUpperBoundScalarField_;
    vtkDoubleArray        **outputProbabilityScalarField_;
    vtkDoubleArray        *outputMeanField_;

    // base code features
    int doIt(vtkDataSet **input,
             vtkDataSet *outputBoundFields,
             vtkDataSet *ouputProbability,
             vtkDataSet *outputMean,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

};

#endif // _TTK_UNCERTAINDATAESTIMATOR_H
