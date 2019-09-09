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

#include <vtkTableWriter.h>
#include <vtkThreshold.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

int main(int argc, char **argv) {

  ttk::CommandLineParser parser;
  ttk::globalDebugLevel_ = 3;

  std::string inputFilePath;

  parser.setArgument("i", &inputFilePath, "Path to input VTU file");
  parser.parse(argc, argv);

  // 1. loading the input data
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader
    = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(inputFilePath.data());

  // 2. computing the persistence curve
  vtkSmartPointer<ttkPersistenceCurve> curve
    = vtkSmartPointer<ttkPersistenceCurve>::New();
  curve->SetInputConnection(reader->GetOutputPort());

  // 3. computing the persitence diagram
  vtkSmartPointer<ttkPersistenceDiagram> diagram
    = vtkSmartPointer<ttkPersistenceDiagram>::New();
  diagram->SetInputConnection(reader->GetOutputPort());

  // 4. selecting the critical point pairs
  vtkSmartPointer<vtkThreshold> criticalPairs
    = vtkSmartPointer<vtkThreshold>::New();
  criticalPairs->SetInputConnection(diagram->GetOutputPort());
  criticalPairs->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "PairIdentifier");
  criticalPairs->ThresholdBetween(-0.1, 999999);

  // 5. selecting the most persistent pairs
  vtkSmartPointer<vtkThreshold> persistentPairs
    = vtkSmartPointer<vtkThreshold>::New();
  persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
  persistentPairs->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
  persistentPairs->ThresholdBetween(0.05, 999999);

  // 6. simplifying the input data to remove non-persistent pairs
  vtkSmartPointer<ttkTopologicalSimplification> topologicalSimplification
    = vtkSmartPointer<ttkTopologicalSimplification>::New();
  topologicalSimplification->SetInputConnection(0, reader->GetOutputPort());
  topologicalSimplification->SetInputConnection(
    1, persistentPairs->GetOutputPort());

  // 7. computing the Morse-Smale complex
  vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex
    = vtkSmartPointer<ttkMorseSmaleComplex>::New();
  morseSmaleComplex->SetInputConnection(
    topologicalSimplification->GetOutputPort());

  // 8. saving the output data
  vtkSmartPointer<vtkTableWriter> curveWriter
    = vtkSmartPointer<vtkTableWriter>::New();
  curveWriter->SetInputConnection(curve->GetOutputPort());
  curveWriter->SetFileName("curve.vtk");
  curveWriter->Write();

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> sepWriter
    = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  sepWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
  sepWriter->SetFileName("separatrices.vtu");
  sepWriter->Write();

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> segWriter
    = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  segWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
  segWriter->SetFileName("segmentation.vtu");
  segWriter->Write();

  return 0;
}

/// @}
