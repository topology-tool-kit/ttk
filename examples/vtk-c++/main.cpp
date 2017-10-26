#include <vtkMorseSmaleComplex.h>
#include <vtkPersistenceCurve.h> 
#include <vtkPersistenceDiagram.h>
#include <vtkTableWriter.h>
#include <vtkThreshold.h>
#include <vtkTopologicalSimplification.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

int main(int argc, char **argv){
  
  // 1. loading the input data
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = 
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName("inputData.vtu");
 
  // 2. computing the persistence curve
  vtkSmartPointer<vtkPersistenceCurve> curve = 
    vtkSmartPointer<vtkPersistenceCurve>::New();
  curve->SetInputConnection(reader->GetOutputPort());
  
  // 3. computing the persitence diagram
  vtkSmartPointer<vtkPersistenceDiagram> diagram = 
    vtkSmartPointer<vtkPersistenceDiagram>::New();
  diagram->SetInputConnection(reader->GetOutputPort());
  
  // 4. selecting the critical point pairs
  vtkSmartPointer<vtkThreshold> criticalPairs = 
    vtkSmartPointer<vtkThreshold>::New();
  criticalPairs->SetInputConnection(diagram->GetOutputPort());
  criticalPairs->SetInputArrayToProcess(0, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_CELLS, "PairIdentifier");
  criticalPairs->ThresholdBetween(-0.1, 999999);
  
  // 5. selecting the most persistent pairs
  vtkSmartPointer<vtkThreshold> persistentPairs = 
    vtkSmartPointer<vtkThreshold>::New();
  persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
  persistentPairs->SetInputArrayToProcess(0, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
  persistentPairs->ThresholdBetween(1.0, 999999);
  
  // 6. simplifying the input data to remove non-persistent pairs
  vtkSmartPointer<vtkTopologicalSimplification> topologicalSimplification =
    vtkSmartPointer<vtkTopologicalSimplification>::New();
  topologicalSimplification->SetInputConnection(0, reader->GetOutputPort());
  topologicalSimplification->SetInputConnection(1, 
    persistentPairs->GetOutputPort());
  
  // 7. computing the Morse-Smale complex
  vtkSmartPointer<vtkMorseSmaleComplex> morseSmaleComplex = 
    vtkSmartPointer<vtkMorseSmaleComplex>::New();
  morseSmaleComplex->SetInputConnection(
    topologicalSimplification->GetOutputPort());
  morseSmaleComplex->SetUseInputOffsetScalarField(true);
  
  // 8. saving the output data
  vtkSmartPointer<vtkTableWriter> curveWriter = 
    vtkSmartPointer<vtkTableWriter>::New();
  curveWriter->SetInputConnection(curve->GetOutputPort());
  curveWriter->SetFileName("curve.vtk");
  curveWriter->Write();
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> sepWriter = 
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  sepWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
  sepWriter->SetFileName("separatrices.vtu");
  sepWriter->Write();
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> segWriter = 
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  segWriter->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
  segWriter->SetFileName("segmentation.vtu");
  segWriter->Write();
  
  return 0;
}