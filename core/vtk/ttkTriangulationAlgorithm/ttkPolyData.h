#ifndef _TTK_POLYDATA_H
#define _TTK_POLYDATA_H

#include <macro.h>
#include <ttkTriangulation.h>

#include <vtkPolyData.h>

#include <ttkTriangulationAlgorithmModule.h>

class TTKTRIANGULATIONALGORITHM_EXPORT ttkPolyData : public ttkTriangulation,
                                                     public vtkPolyData {

public:
  static ttkPolyData *New();
  ttkTypeMacro(ttkPolyData, vtkPolyData);

  void CopyStructure(vtkDataSet *other) override;

  void DeepCopy(vtkDataObject *other) override;

  void ShallowCopy(vtkDataObject *other) override;

protected:
  ttkPolyData();

  ~ttkPolyData() override;
};

#endif /* end of include guard: _TTK_POLYDATA_H */
