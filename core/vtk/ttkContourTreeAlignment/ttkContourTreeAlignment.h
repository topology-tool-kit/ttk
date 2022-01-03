/// \ingroup vtk
/// \class ttkContourTreeAlignment
/// \author Florian Wetzels (f_wetzels13@cs.uni-kl.de), Jonas Lukasczyk
/// (jl@jluk.de) \date 28.01.2020
///
/// \brief TTK VTK-filter that computes an alignment for a multiblock of
/// contourtrees
///
/// VTK wrapping code for the ttk::ContourTreeAlignment package.
///
/// This filter takes a multiblock of unstructured grids, where each block
/// represents a contour tree, and computed the alignment of these contour
/// trees. For each tree, a point array for the scalar vlues, a cell array for
/// the region sizes of arcs and a cell array "for the segmentation ids of arcs
/// are required. Contour trees computed by the FTMTree module fulfill these
/// requirements and are the recommended way to use this filter.
///
/// \param Input vtkMultiBlockDataSet representing the contour trees to align
/// \param Output vtkUnstructuredGrid representing the alignment tree
///
/// The input data array for the scalar values needs to be specified via the
/// standard VTK call vtkAlgorithm::SetInputArrayToProcess() with the following
/// parameters: \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The input data array for the region sizes needs to be specified via the
/// standard VTK call vtkAlgorithm::SetInputArrayToProcess() with the following
/// parameters: \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)R
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 1 (FIXED: cell data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The input data array for the segmentation ids needs to be specified via the
/// standard VTK call vtkAlgorithm::SetInputArrayToProcess() with the following
/// parameters: \param idx 2 (FIXED: the third array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 1 (FIXED: cell data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// \b Related \b publication: \n
/// 'Fuzzy contour trees: Alignment and joint layout of multiple contour trees'
/// Anna Pia Lohfink, Florian Wetzels, Jonas Lukasczyk, Gunther H. Weber, and
/// Christoph Garth. Comput. Graph. Forum, 39(3):343-355, 2020.
///
/// \sa ttk::ContourTreeAlignment
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">
///   Contour Tree Alignment example</a> \n

#pragma once

// VTK Module
#include <ttkContourTreeAlignmentModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <ContourTreeAlignment.h>

class TTKCONTOURTREEALIGNMENT_EXPORT ttkContourTreeAlignment
  : public ttkAlgorithm,
    protected ttk::ContourTreeAlignment {

public:
  // VTK stuff
  static ttkContourTreeAlignment *New();
  vtkTypeMacro(ttkContourTreeAlignment, ttkAlgorithm);

  // auto generated setters and getters
  vtkSetMacro(RandomSeed, int);
  vtkGetMacro(RandomSeed, int);
  vtkSetMacro(MatchTime, bool);
  vtkGetMacro(MatchTime, bool);
  vtkSetMacro(ArcMatchMode, int);
  vtkGetMacro(ArcMatchMode, int);
  vtkSetMacro(AlignmenttreeType, int);
  vtkGetMacro(AlignmenttreeType, int);
  vtkSetMacro(WeightCombinatorialMatch, float);
  vtkGetMacro(WeightCombinatorialMatch, float);
  vtkSetMacro(WeightArcMatch, float);
  vtkGetMacro(WeightArcMatch, float);
  vtkSetMacro(WeightScalarValueMatch, float);
  vtkGetMacro(WeightScalarValueMatch, float);
  vtkSetMacro(ExportJSON, bool);
  vtkGetMacro(ExportJSON, bool);
  vtkSetMacro(ExportPath, const std::string &);
  vtkGetMacro(ExportPath, std::string);

protected:
  // Specify the input data type of each input port
  int FillInputPortInformation(int port, vtkInformation *info) override;

  // Specify the data object type of each output port
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  // Pass VTK data to the base code and convert base code output to VTK
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  // filter constructor and destructor
  ttkContourTreeAlignment();
  ~ttkContourTreeAlignment() override{};

private:
  // filter parameters
  int RandomSeed{1};
  bool MatchTime{false};
  int AlignmenttreeType{0};
  bool ExportJSON{false};
  std::string ExportPath{""};
  int ArcMatchMode{2};
  float WeightCombinatorialMatch{0};
  float WeightArcMatch{1};
  float WeightScalarValueMatch{0};
};
