/// \ingroup vtk
/// \class ttkTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief ttkTriangulation is a wrapper for the ttk::Triangulation class, which
/// provides time and memory efficient traversal methods on triangulations of
/// piecewise linear manifolds. It provides the following features:
///   -# Given a vertex, it provides: the list of edges that are connected to
/// it, the list of its neighbors, its link, its star, etc.
///   -# Given an edge, it provides: its vertices, its star, etc.
///   -# Given a triangle, its provides: its vertices, its edges, etc.
///   -# Given a tetrahedron, its provides: its vertices, its edges, its
/// neighbor tetrahedra, etc.
///   -# Given a triangulation, it provides: its list of vertices, edges,
/// triangles and tetrahedra.
///
/// ttk::Triangulation implements faster accesses than, more general-purpose,
/// competing VTK data structures such as the vtkUnstructuredGrid class.
///
/// ttk::Triangulation supports both explicit and implicit triangulations:
///   -# Explicit triangulations with vtkUnstructuredGrid or vtkPolyData
/// objects: Given a vtkUnstructuredGrid or a vtkPolyData representing a valid
/// triangulation, ttk:Triangulation provides time efficient accesses
/// (requiring adequate pre-processing, see the ttk::Triangulation class
/// documentation).
///   -# Implicit triangulations with vtkImageData objects: Given a vtkImageData
/// representing a regular grid, ttk::Triangulation will perform an implicit
/// triangulation of the grid, enabling both time and memory efficient
/// traversals of triangulations of regular grids.
///
/// Apart from pre-processes, ttk::Triangulation requires no memory overhead in
/// addition to the input vtkUnstructuredGrid, vtkPolyData or vtkImageData
/// objects (the actual data is passed through pointers). This also means that
/// the input VTK objects (vtkUnstructuredGrid, vtkPolyData or vtkImageData)
/// must persist as long as their corresponding ttk::Triangulation object
/// (unspecified behavior otherwise).
///
/// \note
/// Only pre-process the information you need! See the
/// ttk::Triangulation class documentation.
/// \sa ttk::Triangulation

#ifndef _TTK_TRIANGULATION_H
#define _TTK_TRIANGULATION_H

// c++ includes
#include <algorithm>

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

// VTK includes
#include <vtkCellArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#ifndef TTK_PLUGIN
class VTKCOMMONDATAMODEL_EXPORT ttkTriangulation : public ttk::Debug {
#else
class ttkTriangulation : public ttk::Debug {
#endif

public:
  ttkTriangulation();
  ~ttkTriangulation();

  /// Allocates the memory for the internal ttk::Triangulation data-structure.
  /// \warning Internal TTK usage only.
  int allocate();

  /// Retrieves a pointer to the internal ttk::Triangulation of the current
  /// object. Traversal can only be performed on a ttk::Triangulation object.
  /// \return  Returns a pointer to a valid ttk::Triangulation object upon
  /// success, NULL otherwise.
  /// \sa ttk::Triangulation
  ttk::Triangulation *getTriangulation() {
    return triangulation_;
  };

  /// Retrieves a pointer to a ttk::Triangulation object from a vtkDataSet.
  /// \warning This function will return a non null pointer if and only if
  /// the VTK data-set has undergone a ttkTriangulationFilter (i.e. if it
  /// is the direct output of a ttkTriangulationFilter or if a
  /// ttkTriangulationFilter has been called in an earlier stage of the
  /// pipeline).
  /// \param dataSet Input VTK data-set.
  /// \return Returns a pointer to a valid ttk::Triangulation object upon
  /// success, NULL otherwise.
  /// \sa ttk::Triangulation
  /// \sa ttkWrapper
  /// \sa ttkTriangulationFilter
  static ttk::Triangulation *getTriangulation(vtkDataSet *dataSet);

  /// Translates the current triangulation into a vtkUnstructuredGrid object.
  /// This function, used in conjunction with a vtkXMLUnstructuredGridWriter
  /// object, can be useful to store the current triangulation on disk as a
  /// VTU file.
  /// \return Returns a valid pointer to a vtkUnstructuredGrid object
  /// representing the current triangulation.
  ///
  /// \warning If the encapsulated ttk::Triangulation object is representing
  /// implicitly a regular grid, this grid will be explicitly converted into
  /// a triangulation in the output vtkUnstructuredGrid object. This is
  /// likely to be very memory-expensive.
  vtkUnstructuredGrid *getVtkUnstructuredGrid();

  /// Check if the connectivity of the input data-set \p dataSet changed
  /// since the last time the VTK object \p callingObject has been modified.
  /// This function is particularly useful in vtkWrappers to check if the
  /// input geometry changed (requiring a reset of the computation or not).
  /// \param dataSet Input VTK data-set.
  /// \param callingObject Calling VTK object (typically a vtkWrapper)
  /// \return Returns 0 upon success, negative values otherwise.
  /// \sa ttkFTMTree
  /// \sa vtkMorseSmaleComplex
  /// \sa vtkReebSpace
  static bool hasChangedConnectivity(ttk::Triangulation *triangulation,
                                     vtkDataSet *dataSet,
                                     vtkObject *callingObject);

  /// Specify the input VTK object representing a triangulation or a regular
  /// grid.
  /// \param dataSet Input VTK object (vtkImageData, vtkPolyData or
  /// vtkUnstructuredGrid).
  /// \return Returns 0 upon success, negative values otherwise.
  ///
  /// \warning If the internal ttk::Triangulation object is already
  /// representing a valid triangulation, this information will be
  /// over-written (which means that pre-processing functions should
  /// be called again).
  ///
  /// \warning For memory-efficiency purposes, this function makes NO DATA
  /// COPY. This means that, if the object \p dataSet is destroyed after this
  /// function, the behavior of the internal ttk::Triangulation will be
  /// unspecified (well, this is a nice way to say it's gonna crash).
  int setInputData(vtkDataSet *dataSet);

protected:
  int deepCopy(vtkDataObject *other);

  int shallowCopy(vtkDataObject *other);

  bool hasAllocated_;

  vtkDataSet *inputDataSet_;
  vtkSmartPointer<vtkPoints> vtkPoints_;
  vtkSmartPointer<vtkUnstructuredGrid> vtkUnstructuredGrid_;

  ttk::Triangulation *triangulation_;

private:
};

// Internal things to allow the ttkTriangulation to travel through a VTK
// pipeline.
#define ttkTypeMacro(thisClass, superClass)                      \
protected:                                                       \
  const char *GetClassNameInternal() const {                     \
    return #superClass;                                          \
  }                                                              \
                                                                 \
public:                                                          \
  typedef superClass Superclass;                                 \
  static bool IsTypeOf(const char *type) {                       \
    if(!strcmp("superClass", type)) {                            \
      return 1;                                                  \
    }                                                            \
    return superClass::IsTypeOf(type);                           \
  }                                                              \
  int IsA(const char *type) {                                    \
    return this->thisClass::IsTypeOf(type);                      \
  }                                                              \
  static thisClass *SafeDownCast(vtkObjectBase *o) {             \
    if((o) && (o->IsA("thisClass"))) {                           \
      return static_cast<thisClass *>(o);                        \
    }                                                            \
    return NULL;                                                 \
  }                                                              \
  thisClass *NewInstance() const {                               \
    return thisClass::SafeDownCast(this->NewInstanceInternal()); \
  }                                                              \
                                                                 \
protected:                                                       \
  vtkObjectBase *NewInstanceInternal() const {                   \
    return thisClass::New();                                     \
  }                                                              \
                                                                 \
public:

#ifndef TTK_PLUGIN
class VTKCOMMONDATAMODEL_EXPORT ttkUnstructuredGrid :
#else
class ttkUnstructuredGrid :
#endif
  public ttkTriangulation,
  public vtkUnstructuredGrid {

public:
  static ttkUnstructuredGrid *New();
  ttkTypeMacro(ttkUnstructuredGrid, vtkUnstructuredGrid);

  void CopyStructure(vtkDataSet *other);

  void DeepCopy(vtkDataObject *other);

  void ShallowCopy(vtkDataObject *other);

protected:
  ttkUnstructuredGrid();

  ~ttkUnstructuredGrid();
};

#ifndef TTK_PLUGIN
class VTKCOMMONDATAMODEL_EXPORT ttkImageData :
#else
class ttkImageData :
#endif
  public ttkTriangulation,
  public vtkImageData {

public:
  static ttkImageData *New();
  ttkTypeMacro(ttkImageData, vtkImageData);

  void CopyStructure(vtkDataSet *other);

  void DeepCopy(vtkDataObject *other);

  void ShallowCopy(vtkDataObject *other);

protected:
  ttkImageData();

  ~ttkImageData();
};

#ifndef TTK_PLUGIN
class VTKCOMMONDATAMODEL_EXPORT ttkPolyData :
#else
class ttkPolyData :
#endif
  public ttkTriangulation,
  public vtkPolyData {

public:
  static ttkPolyData *New();
  ttkTypeMacro(ttkPolyData, vtkPolyData);

  void CopyStructure(vtkDataSet *other);

  void DeepCopy(vtkDataObject *other);

  void ShallowCopy(vtkDataObject *other);

protected:
  ttkPolyData();

  ~ttkPolyData();
};

#endif // _TTK_TRIANGULATION_H
