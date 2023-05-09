/// \ingroup vtk
/// \class ttkPathCompression
/// \author Robin Maack <maack@rptu.de>
/// \date January 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::PathCompression module.
///
/// This VTK filter uses the ttk::PathCompression module to compute a
/// Morse-Smale segmentation using path compression
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/PathCompression/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
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

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

protected:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputScalars,
               const SimplexId *const inputOffsets,
               const triangulationType &triangulation);

  ttkPathCompression();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputOffsetScalarField{};
  OutputManifold segmentations_{};
};
