#ifndef _TTK_IMAGEDATA_H
#define _TTK_IMAGEDATA_H

#include <macro.h>
#include <ttkTriangulation.h>

#include <vtkImageData.h>

#include <ttkTriangulationAlgorithmModule.h>

class TTKTRIANGULATIONALGORITHM_EXPORT ttkImageData : public ttkTriangulation,
                                                      public vtkImageData {

public:
  static ttkImageData *New();
  ttkTypeMacro(ttkImageData, vtkImageData);

  void CopyStructure(vtkDataSet *other) override;

  void DeepCopy(vtkDataObject *other) override;

  void ShallowCopy(vtkDataObject *other) override;

protected:
  ttkImageData();

  ~ttkImageData() override;
};

#endif /* end of include guard: _TTK_IMAGEDATA_H */
