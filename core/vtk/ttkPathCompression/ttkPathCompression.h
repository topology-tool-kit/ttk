/// \ingroup vtk
/// \class ttkPathCompression
/// \author Robin Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::PathCompression module.
///
/// Given an input order field, this class computes its ascending and descending
/// segmentation by assigning every vertex to its minimum or maximum in gradient
/// or inverse gradient direction. For convienience a hash (no hash collision
/// detection) of both segmentations can be created to represent the Morse-Smale
/// segmentation.
///
/// \param Input vtkDataSet containing the input scalar field as point data
/// \param Output vtkDataSet containing the output segmentations as point data
/// arrays
///
/// This filter can be used like any other VTK filter (for instance, by using
/// the sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// See the corresponding standalone program for a usage example:
///   - standalone/PathCompression/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/TODO/">TODO
///   example</a> \n
///
/// \sa ttk::PathCompression
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkPathCompressionModule.h>

// ttk code includes
#include <PathCompression.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkPolyData;

class TTKPATHCOMPRESSION_EXPORT ttkPathCompression
  : public ttkAlgorithm,
    protected ttk::PathCompression {

public:
  static ttkPathCompression *New();

  vtkTypeMacro(ttkPathCompression, ttkAlgorithm);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeMSSegmentationHash, bool);
  vtkGetMacro(ComputeMSSegmentationHash, bool);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

protected:
  template <typename triangulationType>
  int dispatch(const SimplexId *const inputOffsets,
               const triangulationType &triangulation);

  ttkPathCompression();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputOffsetScalarField{};
  OutputSegmentation segmentations_{};
};
