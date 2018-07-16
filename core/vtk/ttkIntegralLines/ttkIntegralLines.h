/// \ingroup vtk
/// \class ttkIntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK VTK-filter for the computation of edge-based integral lines of 
/// the gradient of an input scalar field. 
/// 
/// The filter takes on its input a scalar field attached as point data to an 
/// input geometry (either 2D or 3D, either regular grids or triangulations) 
/// and computes the forward or backward integral lines along the edges of the
/// input mesh, given a list of input sources. 
/// The sources are specified with a vtkPointSet on which is attached as point 
/// data a scalar field that represent the vertex identifiers of the sources in 
/// the input geometry.
///
/// \param Input0 Input scalar field, either 2D or 3D, either regular grid or 
/// triangulation (vtkDataSet)
/// \param Input1 Input sources (vtkPointSet)
/// \param Output Output integral lines (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::IntegralLines
/// \sa vtkIdentifiers

#ifndef _TTK_DISCRETESTREAMLINE_H
#define _TTK_DISCRETESTREAMLINE_H

// VTK includes
#include<vtkLine.h>
#include<vtkDataArray.h>
#include<vtkDoubleArray.h>
#include<vtkDataSet.h>
#include<vtkPointData.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>

// ttk code includes
#include<IntegralLines.h>
#include<ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkIntegralLines
#else
class ttkIntegralLines
#endif
: public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkIntegralLines* New();

    vtkTypeMacro(ttkIntegralLines, vtkDataSetAlgorithm);

    vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }

    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

    vtkGetMacro(Direction, int);
    vtkSetMacro(Direction, int);

    vtkSetMacro(OutputScalarFieldType, int);
    vtkGetMacro(OutputScalarFieldType, int);

    vtkSetMacro(VertexIdentifierScalarFieldName, std::string);
    vtkGetMacro(VertexIdentifierScalarFieldName, std::string);

    vtkSetMacro(UseOffsetScalarField, int);
    vtkGetMacro(UseOffsetScalarField, int);

    vtkSetMacro(OffsetScalarFieldName, std::string);
    vtkGetMacro(OffsetScalarFieldName, std::string);

    int getTriangulation(vtkDataSet* input);
    int getScalars(vtkDataSet* input);
    int getOffsets(vtkDataSet* input);
    int getIdentifiers(vtkPointSet* input);
    int getTrajectories(vtkDataSet* input, std::vector<std::vector<ttkIdType>>& trajectories, vtkUnstructuredGrid* output);

  protected:

    ttkIntegralLines();
    ~ttkIntegralLines();
    
    TTK_SETUP();

    int FillInputPortInformation(int port, vtkInformation* info);
    int FillOutputPortInformation(int port, vtkInformation* info);

  private:

    bool hasUpdatedMesh_;
    std::string ScalarField;
    int Direction;
    int OutputScalarFieldType;
    std::string VertexIdentifierScalarFieldName;
    int UseOffsetScalarField;
    std::string OffsetScalarFieldName;

    ttk::Triangulation *triangulation_;
    ttk::IntegralLines integralLines_;
    vtkDataArray* inputScalars_;
    vtkDataArray* offsets_;
    vtkDataArray* inputOffsets_;
    vtkDataArray* identifiers_;

};

#endif // _TTK_DISCRETESTREAMLINE_H
