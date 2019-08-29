## TTK - ChangeLog
=
### ttk.git

### 0.9.8
- Updates for ParaView-5.6.1
- ContourAroundPoint (EnvirVis 2019 paper)
- Task-based Parallel Reeb Graphs with Dynamic ST-Trees (EGPGV 2019 paper)
- Nested Tracking Graphs (EuroVis 2017 paper)
- Morse-Smale quadrangulation (SIGGRAPH 2006 paper)
- Harmonic scalar field design (SMI 2009 paper)
- Eigen functions of the Laplace-Beltrami operator (cotan weights)
- Eigen and spectra support
- Barycentric subdivision
- Updated python examples
- Docker support
- Planar graph layout module
- Improved tracking from overlap (modular, templated, streaming support)
- New identifier sorter (for example to select the N most persistent features)
- Memory footprint improvement for the discrete gradient
- New unified plugin library
- Improved memory
- Bug fixes
- Clang-formatting

### 0.9.7
- Updates for ParaView-5.6.0
- Tracking from overlap
- TTK pipeline filters (ForEachRow)
- Depth image based geometry approximation (IEEE LDAV 2018 paper)
- Advanced cinema+sqlite3 support
- Lifted Wasserstein Matcher for Topology Time-Tracking (IEEE LDAV 2018 paper)
- Dimension reduction for high dimensional data (scikit-learn integration)
- Automatic guesses for ttk generated fields
- Improved type consistency
- Improved ParaView GUIs
- DataSet Interpolator (interpolate anything onto anything)
- Improved triangulation request
- Bug fixes

### 0.9.6
- Fixed major performance bug with discrete gradient (thanks to Attila Gyulassy)
- TDA aware compression (with ZLIB and ZFP support)
- New addressing scheme (allowing up to 64 bit ids), for extreme size datasets
- Automatic offset selection
- Updates for ParaView-5.5.2
- Various bug fixes
- Improved examples
- Basic windows continuous integration support
- 1000th commit!

### 0.9.5
- Updates for ParaView-5.5.0
- Bug fixes
- Removed global namespace usage in headers
- Fixes for gcc 4
- Fixes for clang
- Improved triangulation consistency in 2D and 1D
- Improved consistency in morse-smale complex computation
- Manifold check module
- Point merger
- Minimalist OBJ writer
- Support for double and single precision for point coordinates

### 0.9.4
- Windows support
- OFF reader/writer
- Field selector
- Bug fixes

### 0.9.3
- Updates for ParaView-5.4.1
- New cmake structure, find_package(TTKBase) find_package(TTKVTK) now available
- New simple python, vtk-c++ and plain c++ examples (TTK companion paper teaser)
- Greatly accelerated merge tree features using task-based Fibonacci heaps
- Bottleneck and Wasserstein distances between peristence diagrams
- Lp distance between scalar fields
- boundary mask on scalar field critical points
- optional 3D embedding of the persistence diagram
- optional mask array for constrained smoothing
- bug fixes in DiscreteGradient
- bug fixes in TopologicalSimplification
- code cleaning

### 0.9.2
- Updates for ParaView-5.4.0
- New naming convention for vtkWrappers (name conflicts under MacOs)
- TriangulationRequest module (to locally inspect triangulations in ParaView)
- improved time efficiency for link queries on implicit triangulations
- constructors fixes for VTK wrappers
- removed legacy CMake options
- minor bug fixes

### 0.9.1
- latest triangulation API
- improved module handling scripts
- build fixes (in response to early user feedback)
- updates for VTK-7.1.0 and ParaView-5.3.0
- IdentifierRandomizer module (to shuffle segmentation ids)
