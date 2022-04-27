#!/usr/bin/env python

# /// \ingroup examples
# /// \author Lutz Hofmann <lutz.hofmann@iwr.uni-heidelberg.de>
# /// \date August 2019.
# ///
# /// \brief Minimalist python TTK example pipeline, including:
# ///  -# The computation of a persistence curve
# ///  -# The computation of a persistence diagram
# ///  -# The selection of the most persistent pairs of the diagram
# ///  -# The pre-simplification of the data according to this selection
# ///  -# The computation of the Morse-Smale complex on this simplified data
# ///  -# The storage of the output of this pipeline to disk.
# ///
# /// This example reproduces the Figure 1 of the TTK companion paper:
# /// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
# /// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
# /// of IEEE VIS 2017.

import sys

from vtk import (
    vtkDataObject,
    vtkTableWriter,
    vtkThreshold,
    vtkXMLPolyDataWriter,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
)

from topologytoolkit import (
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
curve.SetInputArrayToProcess(0, 0, 0, 0, "data")
curve.SetDebugLevel(3)

# 3. computing the persitence diagram
diagram = ttkPersistenceDiagram()
diagram.SetInputConnection(reader.GetOutputPort())
diagram.SetInputArrayToProcess(0, 0, 0, 0, "data")
diagram.SetDebugLevel(3)

# 4. selecting the critical point pairs
criticalPairs = vtkThreshold()
criticalPairs.SetInputConnection(diagram.GetOutputPort())
criticalPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "PairIdentifier"
)
criticalPairs.ThresholdBetween(-0.1, 999999)

# 5. selecting the most persistent pairs
persistentPairs = vtkThreshold()
persistentPairs.SetInputConnection(criticalPairs.GetOutputPort())
persistentPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "Persistence"
)
persistentPairs.ThresholdBetween(0.05, 999999)

# 6. simplifying the input data to remove non-persistent pairs
topologicalSimplification = ttkTopologicalSimplification()
topologicalSimplification.SetInputConnection(0, reader.GetOutputPort())
topologicalSimplification.SetInputArrayToProcess(0, 0, 0, 0, "data")
topologicalSimplification.SetInputConnection(1, persistentPairs.GetOutputPort())
topologicalSimplification.SetDebugLevel(3)

# 7. computing the Morse-Smale complex
morseSmaleComplex = ttkMorseSmaleComplex()
morseSmaleComplex.SetInputConnection(topologicalSimplification.GetOutputPort())
morseSmaleComplex.SetInputArrayToProcess(0, 0, 0, 0, "data")
morseSmaleComplex.SetDebugLevel(3)

# 8. saving the output data
curveWriter = vtkTableWriter()
curveWriter.SetInputConnection(curve.GetOutputPort())
curveWriter.SetFileName("curve.vtk")
curveWriter.Write()

sepWriter = vtkXMLPolyDataWriter()
sepWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
sepWriter.SetFileName("separatrices.vtp")
sepWriter.Write()

segWriter = vtkXMLUnstructuredGridWriter()
segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(3))
segWriter.SetFileName("segmentation.vtu")
segWriter.Write()
