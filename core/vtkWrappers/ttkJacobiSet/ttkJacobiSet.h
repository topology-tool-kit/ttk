/// \ingroup vtkWrappers
/// \class ttkJacobiSet
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2015.
///
/// \brief TTK VTK-filter that computes the Jacobi set of a bivariate 
/// volumetric data-set.
///
/// Given a bivariate scalar field defined on a PL 3-manifold, this filter
/// produces the list of Jacobi edges (each entry is a pair given by the edge
/// identifier and the Jacobi edge type).
///
/// The input bivariate data must be provided as two independent scalar fields
/// attached as point data to the input geometry.
///
/// \param Input Input bivariate volumetric data (vtkDataSet)
/// \param Output Output Jacobi set (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Jacobi sets of multiple Morse functions" \n
/// Herbert Edelsbrunner, John Harer \n
/// Foundations of Computational Mathematics. Cambridge University Press, 2002.
///
/// \sa ttk::JacobiSet
/// \sa vtkReebSpace

#ifndef _TTK_JACOBISET_H
#define _TTK_JACOBISET_H

// ttk code includes
#include                  <JacobiSet.h>
#include                  <ttkWrapper.h>

// VTK includes -- to adapt
#include                  <vtkCellArray.h>
#include                  <vtkCellData.h>
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnsignedShortArray.h>
#include                  <vtkUnstructuredGrid.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
class VTKFILTERSCORE_EXPORT ttkJacobiSet 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkJacobiSet* New();
    
    vtkTypeMacro(ttkJacobiSet, vtkDataSetAlgorithm);
    
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
    
    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    vtkSetMacro(Ucomponent, string);
    vtkGetMacro(Ucomponent, string);
    
    vtkSetMacro(Vcomponent, string);
    vtkGetMacro(Vcomponent, string);

    vtkSetMacro(UcomponentId, int);
    vtkGetMacro(UcomponentId, int);
    
    vtkSetMacro(VcomponentId, int);
    vtkGetMacro(VcomponentId, int);

    vtkSetMacro(UoffsetId, int);
    vtkGetMacro(UoffsetId, int);
    
    vtkSetMacro(VoffsetId, int);
    vtkGetMacro(VoffsetId, int);

    vtkGetMacro(PredefinedOffset, bool);
    vtkSetMacro(PredefinedOffset, bool);
    
    vtkGetMacro(OffsetFieldU, string);
    vtkSetMacro(OffsetFieldU, string);
    
    vtkGetMacro(OffsetFieldV, string);
    vtkSetMacro(OffsetFieldV, string);
    
    vtkSetMacro(EdgeIds, bool);
    vtkGetMacro(EdgeIds, bool);
    
    vtkSetMacro(VertexScalars, bool);
    vtkGetMacro(VertexScalars, bool);
    
    int FillOutputPortInformation(int port, vtkInformation *info){
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
      return 1;
    }
    
  protected:
    
    ttkJacobiSet();
    
    ~ttkJacobiSet();
    
    TTK_SETUP();
    
    
  private:
    
    bool                  PredefinedOffset;
    bool                  EdgeIds, VertexScalars;
    int                   UcomponentId, VcomponentId, UoffsetId, VoffsetId;
    string                Ucomponent, Vcomponent, 
                          OffsetFieldU, OffsetFieldV;
    vector<pair<int, int> > edgeList_;
    // for each edge, one skeleton of its triangle fan
    vector<vector<pair<int, int> > > edgeFanLinkEdgeList_;
    // for each edge, the one skeleton of its triangle fan
    vector<vector<long long int> >  edgeFans_;
    vector<pair<int, char> >        jacobiSet_;
    vector<int>           sosOffsetsU_, sosOffsetsV_;
   
    template<class dataTypeU, class dataTypeV> int baseCall(
      vtkDataSet *input,
      vtkDataArray *uField, 
      vtkDataArray *vField);
    
};

#endif // _TTK_JACOBISET_H
