/// \ingroup vtk
/// \class ttkPersistenceDiagramsBarycenter
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
/// \sa ttk::PersistenceDiagramsBarycenter
#ifndef _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H
#define _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H

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

// ttk code includes
#include                  <PersistenceDiagramsBarycenter.h>
#include                  <Wrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceDiagramsBarycenter
#else
class ttkPersistenceDiagramsBarycenter
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkPersistenceDiagramsBarycenter* New();

    vtkTypeMacro(ttkPersistenceDiagramsBarycenter, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
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

    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);


  protected:

    ttkPersistenceDiagramsBarycenter();

    ~ttkPersistenceDiagramsBarycenter();

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);


  private:

    bool                  UseAllCores;
    int                   ThreadNumber;
    std::string                ScalarField;

    // base code features
    int doIt(vtkDataSet **input,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

};

#endif // _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H
