import sys
from vtk import (
    vtkDataObject,
    vtkTableWriter,
    vtkThreshold,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
)
from ttk import (
    ttkMorseSmaleComplex,
    ttkPersistenceCurve,
    ttkPersistenceDiagram,
    ttkTopologicalSimplification,
)

if len(sys.argv) == 2:
    inputFilePath = sys.argv[1]
else:
    print("Missing mandatory argument: Path to input VTU file")
    sys.exit() 

# 1. loading the input data
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(inputFilePath)

# 2. computing the persistence curve
curve = ttkPersistenceCurve()
curve.SetInputConnection(reader.GetOutputPort())

# 3. computing the persitence diagram
diagram = ttkPersistenceDiagram()
diagram.SetInputConnection(reader.GetOutputPort())

# 4. selecting the critical point pairs
criticalPairs = vtkThreshold()
criticalPairs.SetInputConnection(diagram.GetOutputPort())
criticalPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, 'PairIdentifier')
criticalPairs.ThresholdBetween(-0.1, 999999)

# 5. selecting the most persistent pairs
persistentPairs = vtkThreshold()
persistentPairs.SetInputConnection(criticalPairs.GetOutputPort())
persistentPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, 'Persistence')
persistentPairs.ThresholdBetween(0.05, 999999)

# 6. simplifying the input data to remove non-persistent pairs
topologicalSimplification = ttkTopologicalSimplification()
topologicalSimplification.SetInputConnection(0, reader.GetOutputPort())
topologicalSimplification.SetInputConnection(1, persistentPairs.GetOutputPort())

# 7. computing the Morse-Smale complex
morseSmaleComplex = ttkMorseSmaleComplex();
morseSmaleComplex.SetInputConnection(topologicalSimplification.GetOutputPort())

# 8. saving the output data
curveWriter = vtkTableWriter()
curveWriter.SetInputConnection(curve.GetOutputPort())
curveWriter.SetFileName("curve.vtk")
curveWriter.Write()

sepWriter = vtkXMLUnstructuredGridWriter()
sepWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
sepWriter.SetFileName("separatrices.vtu")
sepWriter.Write()

segWriter = vtkXMLUnstructuredGridWriter()
segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(3))
segWriter.SetFileName("segmentation.vtu")
segWriter.Write()
