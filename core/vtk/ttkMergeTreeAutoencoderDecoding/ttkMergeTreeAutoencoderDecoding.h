/// \ingroup vtk
/// \class ttkMergeTreeAutoencoderDecoding
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeAutoencoderDecoding
/// module.
///
/// This VTK filter uses the ttk::MergeTreeAutoencoderDecoding module to compute
/// a decoding of merge trees or persistence diagrams given the parameters of a
/// Wasserstein Auto-Encoder.
///
/// \param Input vtkMultiBlockDataSet Origins
/// \param Input vtkMultiBlockDataSet Bases Axes
/// \param Input vtkMultiBlockDataSet Coefficients
/// \param Output vtkMultiBlockDataSet Trees
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication: \n
/// "Wasserstein Auto-Encoders of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont,  Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2023
///
/// \sa ttk::MergeTreeAutoencoderDecoding
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeWAE/">Merge
///   tree Wasserstein Auto-Encoder example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramWAE/">Persistence
///   Diagram Wasserstein Auto-Encoder example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeAutoencoderDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <MergeTreeAutoencoderDecoding.h>

class TTKMERGETREEAUTOENCODERDECODING_EXPORT ttkMergeTreeAutoencoderDecoding
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeAutoencoderDecoding // and we inherit from the base
                                                // class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeAutoencoderDecoding *New();
  vtkTypeMacro(ttkMergeTreeAutoencoderDecoding, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeAutoencoderDecoding();
  ~ttkMergeTreeAutoencoderDecoding() override = default;

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
