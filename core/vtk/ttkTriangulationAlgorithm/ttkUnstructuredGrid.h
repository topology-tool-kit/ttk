#ifndef _TTK_UNSTRUCTURED_GRID_H
#define _TTK_UNSTRUCTURED_GRID_H

#include <ttkTriangulation.h>
#include <macro.h>

#include <vtkUnstructuredGrid.h>

class ttkUnstructuredGrid : public ttkTriangulation,
                            public vtkUnstructuredGrid {

public:
  static ttkUnstructuredGrid *New();
  ttkTypeMacro(ttkUnstructuredGrid, vtkUnstructuredGrid);

  void CopyStructure(vtkDataSet *other) override;

  void DeepCopy(vtkDataObject *other) override;

  void ShallowCopy(vtkDataObject *other) override;

protected:
  ttkUnstructuredGrid();

  ~ttkUnstructuredGrid() override;
};

#endif /* end of include guard: _TTK_UNSTRUCTURED_GRID_H */
