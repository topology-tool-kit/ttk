#pragma once

#include <CinemaImagingEmbree.h>

class vtkMultiBlockDataSet;
class vtkPointSet;

namespace ttk {
  class ttkCinemaImagingEmbree : public CinemaImagingEmbree {
  public:
    ttkCinemaImagingEmbree();
    ~ttkCinemaImagingEmbree();

    int RenderVTKObject(vtkMultiBlockDataSet *outputImages,

                        vtkPointSet *inputObject,
                        vtkPointSet *inputGrid) const;
  };
}; // namespace ttk