/// \ingroup vtkWrappers
/// \class vtkContourForests
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>.
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016
///
///\brief TTK VTK-filter that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
/// This filter takes a scalar field attached as point data to a geometry
/// (either 2D or 3D, either regular grids or triangulations) and computes its
/// contour tree. Several outputs are produced to encode the nodes of the tree
/// (as points in 3D space), the arcs of the tree and the data segmentation.
///
/// \param Input Input scalar field, either 2D or 3D, either regular
/// grid or triangulation (vtkDataSet)
/// \param Output0 Output nodes (vtkPolyData)
/// \param Output1 Output arcs (vtkPolyData)
/// \param Output2 Output data segmentation (vtkUnstructuredGrid)
/// \param Output3 Output persistence diagram (vtkUnstructuredGrid)
/// \param Output4 Output persistence curve (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa ttk::ContourForests

#ifndef _VTK_CONTOURTREE_H
#define _VTK_CONTOURTREE_H

// vtk wrapper includes
#include <ttkWrapper.h>

// VTK includes
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTable.h>

// base code includes
#include <Geometry.h>

#include "ContourForests.h"
#include "ContourForestsTree.h"
#include "DataTypes.h"

class VTKFILTERSCORE_EXPORT vtkContourForests 
  : public vtkDataSetAlgorithm, public Wrapper{
    
    
  public:
    
    static vtkContourForests* New();

    vtkTypeMacro(vtkContourForests, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    vtkSetMacro(FieldId, int);

    void SetThreadNumber(int threadNumber);
    void SetDebugLevel(int d);
    void SetUseAllCores(bool onOff);
    // end of default ttk setters

    vtkGetMacro(scalarField_, string);
    void SetScalarField(string scalarField);

    vtkGetMacro(useInputOffsetScalarField_, int);
    void SetUseInputOffsetScalarField(bool onOff);

    vtkSetMacro(inputOffsetScalarFieldName_, string);
    vtkGetMacro(inputOffsetScalarFieldName_, string);

    vtkSetMacro(InputOffsetFieldId, int);
    vtkGetMacro(InputOffsetFieldId, int);

    void SetTreeType(int tree);

    void ShowMin(bool state);
    void ShowMax(bool state);
    void ShowSaddle1(bool state);
    void ShowSaddle2(bool state);

    void ShowArc(bool state);
    void SetArcResolution(int arcResolution);
    void SetPartitionNumber(int partitionNum);
    void SetLessPartition(bool l);

    void SetSkeletonSmoothing(double skeletonSmooth);

    void SetSimplificationType(int type);

    void SetSimplificationThreshold(double simplificationThreshold);

    
  protected:
    
    vtkContourForests();
    ~vtkContourForests();
    
    // VTK Interface //
    virtual int FillInputPortInformation(int port, vtkInformation* info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    // Base //
    int vtkDataSetToStdVector(vtkDataSet* input);
    bool isCoincident(float p1[], double p2[]);
    bool isCoincident(double p1[], double p2[]);

    // ContourForestsTree //
    void getTree();
    void updateTree();
    NodeType getNodeType(int id);
    NodeType getNodeType(int id, TreeType type, MergeTree* tree);
    void getCriticalPoints();
    void clearTree();

    // Skeleton //
    void getSkeleton();
    void clearSkeleton();
    void getSkeletonNodes();
    void getSkeletonArcs();
    int getSkeletonScalars(const vector<double>& scalars,
        vector<vector<double> >& skeletonScalars) const;

    // Segmentation //
    void getSegmentation(vtkDataSet* input);
    void clearSegmentation();

    int sample(unsigned int samplingLevel);

    int computeBarycenters();
    void computeSkeleton(unsigned int arcRes);
    void smoothSkeleton(unsigned int skeletonSmoothing);
    void smooth(const int idArc, bool order);

    TTK_PIPELINE_REQUEST();
    TTK_OUTPUT_MANAGEMENT();
   
    void SetThreads();
    
    int doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs);
    
    bool needsToAbort();
      
    int updateProgress(const float &progress);
   
    
  private:
    // Base //
    bool UseAllCores;
    int ThreadNumber;
    int FieldId;
    int InputOffsetFieldId;
    bool isLoaded_;
    bool lessPartition_;
    MergeTree* tree_;
    ContourForests* contourTree_;
    vtkPolyData* skeletonNodes_;
    vtkPolyData* skeletonArcs_;
    vtkDataSet* segmentation_;

    // Void //
    vtkUnstructuredGrid* voidUnstructuredGrid_;
    vtkPolyData* voidPolyData_;

    // Configuration //
    string inputOffsetScalarFieldName_;
    bool useInputOffsetScalarField_;
    bool varyingMesh_;
    bool varyingDataValues_;
    TreeType treeType_;
    string scalarField_;
    bool showMin_;
    bool showMax_;
    bool showSaddle1_;
    bool showSaddle2_;
    bool showArc_;
    unsigned int arcResolution_;
    int partitionNum_;
    unsigned int skeletonSmoothing_;
    int simplificationType_;
    double simplificationThreshold_;
    double simplificationThresholdBuffer_;

    // Computation handles
    bool toUpdateVertexSoSoffsets_;
    bool toComputeContourTree_;
    bool toUpdateTree_;
    bool toComputeSkeleton_;
    bool toComputeSegmentation_;

    // Convenient storage //
    vtkDataArray* vtkInputScalars_;
    double deltaScalar_;
    unsigned int numberOfVertices_;
    Triangulation *triangulation_;
    vector<vector<int>>* vertexNeighbors_;
    vector<int>* vertexSoSoffsets_;
    vector<int>* criticalPoints_;
    vector<double>* vertexScalars_;
    vector<vector<double>>* inputScalars_;
    vector<string>* inputScalarsName_;

    // treeType, SuperArc, several vertices list.
    vector<vector<vector<vector<int>>>>*    samples_;
    vector<vector<vector<vector<double>>>>* barycenters_;

};

#endif  // _VTK_CONTOURTREE_H
