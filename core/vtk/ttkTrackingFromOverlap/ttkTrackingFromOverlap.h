/// \ingroup vtk
/// \class ttkTrackingFromOverlap
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that computes the overlap between labeled
/// vtkPointSets.
///
/// VTK wrapping code for the ttk::TrackingFromOverlap package.
///
/// This filter identifies and tracks labled vtkPointSets across time (and
/// optionally levels) based on spatial overlap, where two points overlap iff
/// their corresponding coordinates are equal. This filter can be executed
/// iteratively and can generate nested tracking graphs.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'.
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///
/// \param Input A \b vtkMultiBlockDataSet that holds the labeled \b
/// vtkPointSets and has one of the following forms:\n{t_0,...,t_n} or
/// {l_0:{t_0,...,t_n}, ... , l_m:{t_0,...,t_n}} where \b t_i is the \b
/// vtkPointSet of timestep \b i, and \b l_j is a \b vtkMultiBlockDataSet that
/// holds all timesteps of level \b j. \param Output A \b vtkUnstructuredGrid
/// that represents the (nested) tracking graph embedded in the spatial domain.
///
/// sa ttk::TrackingFromOverlap

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridAlgorithm.h>

// VTK Module
#include <ttkTrackingFromOverlapModule.h>

// TTK includes
#include <TrackingFromOverlap.h>
#include <ttkAlgorithm.h>

class TTKTRACKINGFROMOVERLAP_EXPORT ttkTrackingFromOverlap
  : public ttkAlgorithm,
    protected ttk::TrackingFromOverlap {

public:
  static ttkTrackingFromOverlap *New();
  vtkTypeMacro(ttkTrackingFromOverlap, ttkAlgorithm);

  vtkSetMacro(LabelFieldName, const std::string &);
  vtkGetMacro(LabelFieldName, std::string);

protected:
  ttkTrackingFromOverlap() {
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  int reset();

  int packInputData(vtkDataObject *inputDataObject,
                    vtkMultiBlockDataSet *packedData) const;
  int packStreamedData(vtkMultiBlockDataSet *streamedData,
                       vtkMultiBlockDataSet *packedData) const;
  int checkData(vtkMultiBlockDataSet *data);

  int storeStreamedData(vtkMultiBlockDataSet *data);
  int computeNodes(vtkMultiBlockDataSet *data);
  int computeTrackingGraphs(vtkMultiBlockDataSet *data);
  int computeNestingTrees(vtkMultiBlockDataSet *data);
  int computeBranches();

  int meshNestedTrackingGraph(vtkDataObject *trackingGraph);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
        info->Append(
          vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
        info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        return 0;
    }
    return 1;
  }

private:
  int LabelDataType;
  std::string LabelFieldName{"RegionId"};

  vtkSmartPointer<vtkMultiBlockDataSet> previousIterationData;

  // Containers for nodes and edges
  std::vector<std::vector<Nodes>> levelTimeNodesMap; // N
  std::vector<std::vector<Edges>> levelTimeEdgesTMap; // E_T
  std::vector<std::vector<Edges>> timeLevelEdgesNMap; // E_N
};
