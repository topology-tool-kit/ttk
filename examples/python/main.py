#!/usr/bin/pvpython
from paraview.simple import *

# 1. loading the input data
inputData = XMLUnstructuredGridReader(FileName=['inputData.vtu'])

# 2. computing the persistence curve
persistenceCurve = TTKPersistenceCurve(inputData)

# 3. computing the persitence diagram
persistenceDiagram = TTKPersistenceDiagram(inputData)

# 4. selecting the critical point pairs
criticalPointPairs = Threshold(persistenceDiagram)
criticalPointPairs.Scalars = ['CELLS', 'PairIdentifier']
criticalPointPairs.ThresholdRange = [-0.1, 999999]

# 5. selecting the most persistent pairs
persistentPairs = Threshold(criticalPointPairs)
persistentPairs.Scalars = ['CELLS', 'Persistence']
persistentPairs.ThresholdRange = [1.0, 999999]

# 6. simplifying the input data to remove non-persistent pairs
topologicalSimplification = TTKTopologicalSimplification(
  Domain=inputData, Constraints=persistentPairs)

# 7. computing the Morse-Smale complex
morseSmaleComplex = TTKMorseSmaleComplex(topologicalSimplification)
morseSmaleComplex.UseInputOffsetField = 1

# 8. saving the output data
SaveData('curve.csv', OutputPort(persistenceCurve, 3))
SaveData('separatrices.vtu', CleantoGrid(OutputPort(morseSmaleComplex, 1)))
SaveData('segmentation.vtu', CleantoGrid(OutputPort(morseSmaleComplex, 3)))
