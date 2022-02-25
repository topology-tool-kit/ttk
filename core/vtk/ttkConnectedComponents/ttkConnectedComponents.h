/// \ingroup vtk
/// \class ttkConnectedComponents
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 22.02.2022
///
/// \brief TTK VTK-filter that computes connected components based on a scalar
/// field.
///
/// VTK wrapping code for the ttk::ConnectedComponents package.
///
/// This filter consumes a scalar field with a feature mask and computes for
/// each edge connected group of vertices with a non-background mask value a
/// so-called connected component via flood-filling, where the backgroud is
/// masked with values smaller-equal zero. The computed components store the
/// size, seed, and center of mass of each component. The flag
/// UseSeedIdAsComponentId controls if the resulting segmentation is either
/// labeled by the index of the component, or by its seed location (which can be
/// used as a deterministic component label).
///
/// The input data array that contains the feature mask needs to be specified
/// via the standard VTK call SetInputArrayToProcess() with the following
/// parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the VTK array)
///
/// \sa ttk::ConnectedComponents

#pragma once

// VTK Module
#include <ttkConnectedComponentsModule.h>

// TTK Include
#include <ConnectedComponents.h>
#include <ttkAlgorithm.h>

class TTKCONNECTEDCOMPONENTS_EXPORT ttkConnectedComponents
  : public ttkAlgorithm,
    protected ttk::ConnectedComponents {

private:
  double BackgroundThreshold{0.0};
  bool AugmentSegmentationWithComponentSize{false};

public:
  vtkSetMacro(BackgroundThreshold, double);
  vtkGetMacro(BackgroundThreshold, double);
  vtkSetMacro(AugmentSegmentationWithComponentSize, bool);
  vtkGetMacro(AugmentSegmentationWithComponentSize, bool);

  static ttkConnectedComponents *New();
  vtkTypeMacro(ttkConnectedComponents, ttkAlgorithm);

protected:
  ttkConnectedComponents();
  virtual ~ttkConnectedComponents() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
