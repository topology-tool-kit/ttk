#!/usr/bin/env pvpython

#/// \ingroup examples
#/// \author Julien Tierny <julien.tierny@lip6.fr>
#/// \date October 2017.
#/// 
#/// \brief Minimalist python TTK example pipeline, including:
#///  -# The computation of a persistence curve
#///  -# The computation of a persistence diagram
#///  -# The selection of the most persistent pairs of the diagram
#///  -# The pre-simplification of the data according to this selection
#///  -# The computation of the Morse-Smale complex on this simplified data
#///  -# The storage of the output of this pipeline to disk.
#///
#/// This example reproduces the Figure 1 of the TTK companion paper:
#/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
#/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
#/// of IEEE VIS 2017.

from paraview.simple import *

if len(sys.argv) == 2:
    inputFilePath = sys.argv[1]
else:
    print("Missing mandatory argument: Path to input VTU file")
    sys.exit() 

# 1. loading the input data
inputData = XMLUnstructuredGridReader(FileName=[inputFilePath])

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
persistentPairs.ThresholdRange = [0.05, 999999]

# 6. simplifying the input data to remove non-persistent pairs
topologicalSimplification = TTKTopologicalSimplification(
  Domain=inputData, Constraints=persistentPairs)

# 7. computing the Morse-Smale complex
morseSmaleComplex = TTKMorseSmaleComplex(topologicalSimplification)

# 8. saving the output data
SaveData('curve.vtk', OutputPort(persistenceCurve, 3))
SaveData('separatrices.vtu', CleantoGrid(OutputPort(morseSmaleComplex, 1)))
SaveData('segmentation.vtu', CleantoGrid(OutputPort(morseSmaleComplex, 3)))
