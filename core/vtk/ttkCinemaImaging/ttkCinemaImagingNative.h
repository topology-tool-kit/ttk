#pragma once

#include <CinemaImagingNative.h>

class vtkMultiBlockDataSet;
class vtkPointSet;

namespace ttk {
  class ttkCinemaImagingNative : public CinemaImagingNative {
  public:
    ttkCinemaImagingNative();
    ~ttkCinemaImagingNative() override;

    int RenderVTKObject(vtkMultiBlockDataSet *outputImages,

                        vtkPointSet *inputObject,
                        vtkPointSet *inputGrid) const;
  };
}; // namespace ttk