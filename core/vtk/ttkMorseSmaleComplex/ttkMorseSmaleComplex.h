/// \ingroup vtk
/// \class ttkMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK VTK-filter that wraps the morseSmaleComplex processing package.
///
/// TTK module for the computation of Morse-Smale complexes.
/// Morse-Smale complexes are useful topological abstractions of scalar 
/// fields for data segmentation, feature extraction, etc.
/// 
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \param Input Input scalar field, defined as a point data scalar field 
/// attached to a geometry, either 2D or 3D, either regular grid or 
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkUnstructuredGrid)
/// \param Output1 Output 1-separatrices (vtkUnstructuredGrid)
/// \param Output2 Output 2-separatrices (vtkUnstructuredGrid)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::MorseSmaleComplex
///
#ifndef _TTK_MORSESMALECOMPLEX_H
#define _TTK_MORSESMALECOMPLEX_H

#ifndef _MSC_VER
// ttk code includes
#include<MorseSmaleComplex.h>
#include<ttkWrapper.h>

// VTK includes -- to adapt
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkInformationVector.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkCellData.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#else
// VTK includes -- to adapt
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkInformationVector.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkCellData.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>

// ttk code includes
#include<MorseSmaleComplex.h>
#include<ttkWrapper.h>
#endif

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMorseSmaleComplex
#else
class ttkMorseSmaleComplex
#endif
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkMorseSmaleComplex* New();

    vtkTypeMacro(ttkMorseSmaleComplex, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    vtkSetMacro(ScalarFieldId, int);
    vtkGetMacro(ScalarFieldId, int);

    vtkSetMacro(OffsetFieldId, int);
    vtkGetMacro(OffsetFieldId, int);

    vtkSetMacro(UseInputOffsetScalarField, int);
    vtkGetMacro(UseInputOffsetScalarField, int);

    vtkSetMacro(InputOffsetScalarFieldName, string);
    vtkGetMacro(InputOffsetScalarFieldName, string);

    vtkSetMacro(IterationThreshold, int);
    vtkGetMacro(IterationThreshold, int);

    vtkSetMacro(ReverseSaddleMaximumConnection, int);
    vtkGetMacro(ReverseSaddleMaximumConnection, int);

    vtkSetMacro(ReverseSaddleSaddleConnection, int);
    vtkGetMacro(ReverseSaddleSaddleConnection, int);

    vtkSetMacro(ComputeAscendingSeparatrices1, int);
    vtkGetMacro(ComputeAscendingSeparatrices1, int);

    vtkSetMacro(ComputeDescendingSeparatrices1, int);
    vtkGetMacro(ComputeDescendingSeparatrices1, int);

    vtkSetMacro(ComputeSaddleConnectors, int);
    vtkGetMacro(ComputeSaddleConnectors, int);

    vtkSetMacro(ComputeAscendingSeparatrices2, int);
    vtkGetMacro(ComputeAscendingSeparatrices2, int);

    vtkSetMacro(ComputeDescendingSeparatrices2, int);
    vtkGetMacro(ComputeDescendingSeparatrices2, int);

    vtkSetMacro(ComputeAscendingSegmentation, int);
    vtkGetMacro(ComputeAscendingSegmentation, int);

    vtkSetMacro(ComputeDescendingSegmentation, int);
    vtkGetMacro(ComputeDescendingSegmentation, int);

    vtkSetMacro(ComputeFinalSegmentation, int);
    vtkGetMacro(ComputeFinalSegmentation, int);

    vtkSetMacro(ReturnSaddleConnectors, int);
    vtkGetMacro(ReturnSaddleConnectors, int);

    vtkSetMacro(SaddleConnectorsPersistenceThreshold, double);
    vtkGetMacro(SaddleConnectorsPersistenceThreshold, double);

    int setupTriangulation(vtkDataSet* input);
    vtkDataArray* getScalars(vtkDataSet* input);
    vtkDataArray* getOffsets(vtkDataSet* input);

  protected:

    ttkMorseSmaleComplex();
    ~ttkMorseSmaleComplex();

    TTK_SETUP();

    virtual int FillInputPortInformation(int port, vtkInformation* info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

  private:
    string ScalarField;
    string InputOffsetScalarFieldName;
    bool UseInputOffsetScalarField;
    int IterationThreshold;
    bool ReverseSaddleMaximumConnection;
    bool ReverseSaddleSaddleConnection;
    bool ComputeAscendingSeparatrices1;
    bool ComputeDescendingSeparatrices1;
    bool ComputeSaddleConnectors;
    bool ComputeAscendingSeparatrices2;
    bool ComputeDescendingSeparatrices2;
    bool ComputeAscendingSegmentation;
    bool ComputeDescendingSegmentation;
    bool ComputeFinalSegmentation;
    int ScalarFieldId;
    int OffsetFieldId;
    int ReturnSaddleConnectors;
    double SaddleConnectorsPersistenceThreshold;

    MorseSmaleComplex morseSmaleComplex_;
    Triangulation *triangulation_;
    vtkIntArray* defaultOffsets_;
    bool hasUpdatedMesh_;
};

#endif // _TTK_MORSESMALECOMPLEX_H
