/// \defgroup examples examples
/// \brief The Topology ToolKit - Example programs.
/// @{
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2017.
///
/// \brief Minimalist VTK-based TTK example pipeline, including:
///  -# The computation of a persistence curve
///  -# The computation of a persistence diagram
///  -# The selection of the most persistent pairs of the diagram
///  -# The pre-simplification of the data according to this selection
///  -# The computation of the Morse-Smale complex on this simplified data
///  -# The storage of the output of this pipeline to disk.
///
/// This example reproduces the Figure 1 of the TTK companion paper:
/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
/// of IEEE VIS 2017.

#include <CommandLineParser.h>

#include <ttkMorseSmaleComplex.h>
#include <ttkPersistenceCurve.h>
#include <ttkPersistenceDiagram.h>
#include <ttkTopologicalSimplification.h>

#include <vtkNew.h>
#include <vtkTableWriter.h>
#include <vtkThreshold.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

int main(int argc, char **argv) {

  ttk::CommandLineParser parser;
  ttk::globalDebugLevel_ = 3;

  std::string inputFilePath;

  parser.setArgument("i", &inputFilePath, "Path to input VTU file");
  parser.parse(argc, argv);

  // 1. loading the input data
  vtkNew<vtkXMLUnstructuredGridReader> reader{};
  reader->SetFileName(inputFilePath.data());

  // 2. computing the persistence curve
  vtkNew<ttkPersistenceCurve> curve{};
  curve->SetInputConnection(reader->GetOutputPort());
  curve->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "data");

  // 3. computing the persitence diagram
  vtkNew<ttkPersistenceDiagram> diagram{};
  diagram->SetInputConnection(reader->GetOutputPort());
  diagram->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "data");

  // 4. selecting the critical point pairs
  vtkNew<vtkThreshold> criticalPairs{};
  criticalPairs->SetInputConnection(diagram->GetOutputPort());
  criticalPairs->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "PairIdentifier");
  criticalPairs->ThresholdBetween(-0.1, 999999);

  // 5. selecting the most persistent pairs
  vtkNew<vtkThreshold> persistentPairs{};
  persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
  persistentPairs->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
  persistentPairs->ThresholdBetween(0.05, 999999);

  // 6. simplifying the input data to remove non-persistent pairs
  vtkNew<ttkTopologicalSimplification> topologicalSimplification{};
  topologicalSimplification->SetInputConnection(0, reader->GetOutputPort());
  topologicalSimplification->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "data");
  topologicalSimplification->SetInputConnection(
    1, persistentPairs->GetOutputPort());

  // 7. computing the Morse-Smale complex
  vtkNew<ttkMorseSmaleComplex> morseSmaleComplex{};
  morseSmaleComplex->SetInputConnection(
    topologicalSimplification->GetOutputPort());
  morseSmaleComplex->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "data");

  // 8. saving the output data
  vtkNew<vtkTableWriter> curveWriter{};
  curveWriter->SetInputConnection(curve->GetOutputPort());
  curveWriter->SetFileName("curve.vtk");
  curveWriter->Write();

  vtkNew<vtkXMLPolyDataWriter> sepWriter{};
  sepWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
  sepWriter->SetFileName("separatrices.vtp");
  sepWriter->Write();

  vtkNew<vtkXMLUnstructuredGridWriter> segWriter{};
  segWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
  segWriter->SetFileName("segmentation.vtu");
  segWriter->Write();

  return 0;
}

/// @}
