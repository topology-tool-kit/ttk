/// \ingroup vtk
/// \class ttkSignedDistanceField
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Mohamed Amine Kissi <mohamed.kissi@lip6.fr>
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date August 2023
///
/// \brief This filter computes a signed distance field given a surface in
/// input.
///
/// It implements three backends (accelerated with a BVH data structure):
///
/// - The exact backend
///
/// - The fast marching backend, this method simulates a "wave" that move
/// starting from the input surface and solve the eikonal equation vertex by
/// vertex to approximate the signed distance field. It corresponds to the the
/// following reference:
///
/// J.A. Sethian.
/// A Fast Marching Level Set Method for Monotonically Advancing Fronts,
/// Proc. Natl. Acad. Sci., 93, 4, pp.1591--1595, 1996
///
/// - The fast marching band backend, a variant of the fast marching for which
/// all the vertices being not yet updated and nearest the input surface are
/// updated (instead of just one in the fast marching backend). It results in a
/// faster method (due to parallelism) but is a rougher approximation.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topologicalOptimization_pegasus/">Topological
///   Optimization for Pegasus Genus Repair example</a>\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topologicalOptimization_torus/">Topological
///   Optimization for Torus Repair example</a>\n

#pragma once

// VTK Module
#include <ttkSignedDistanceFieldModule.h>

// ttk code includes
#include <SignedDistanceField.h>
#include <ttkAlgorithm.h>

class vtkImageData;

class TTKSIGNEDDISTANCEFIELD_EXPORT ttkSignedDistanceField
  : public ttkAlgorithm,
    protected ttk::SignedDistanceField {
public:
  static ttkSignedDistanceField *New();
  vtkTypeMacro(ttkSignedDistanceField, ttkAlgorithm);

  ///@{
  /**
   * Set/Get sampling dimension along each axis. Default will be [10,10,10]
   */
  vtkSetVector3Macro(SamplingDimensions, int);
  vtkGetVector3Macro(SamplingDimensions, int);
  ///@}

  vtkSetMacro(ExpandBox, bool);
  vtkGetMacro(ExpandBox, bool);

  vtkSetMacro(Backend, int);
  vtkGetMacro(Backend, int);

  vtkSetMacro(FastMarchingOrder, int);
  vtkGetMacro(FastMarchingOrder, int);

  /**
   * Get the output data for this algorithm.
   */
  vtkImageData *GetOutput();

protected:
  ttkSignedDistanceField();

  // Usual data generation method
  vtkTypeBool ProcessRequest(vtkInformation *,
                             vtkInformationVector **,
                             vtkInformationVector *) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  int RequestInformation(vtkInformation *,
                         vtkInformationVector **,
                         vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;
  int FillInputPortInformation(int, vtkInformation *) override;
  int FillOutputPortInformation(int, vtkInformation *) override;

  void computeOutputInformation(vtkInformationVector **inputVector);

  int SamplingDimensions[3];
  bool ExpandBox = true;
  int Backend = 0;
  int FastMarchingOrder = 1;

private:
  std::array<int, 6> DataExtent{0, 0, 0, 0, 0, 0};
  std::array<double, 3> Origin{0.0, 0.0, 0.0};
};
