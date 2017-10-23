# libraries

- [x] make sure only the paraview plugins build against the VTK shipped with PV and not everything.
- [x] rename installation dir to make them consistent with the code source tree (or the other way around)

# CMake TTK messages
- [ ] put the cmake ttk messages back in (useful for debug of low level stuff)

# CMake bundles
- [ ] automatic creation based on the listing of the vtkWrappers directory

# Accessible CMake features
- [x] make ParaView_DIR accessible through cmake config (appears when paraview enabled)
- [x] make VTK_DIR accessible through cmake config
- [x] make TTK_PLUGIN_INSTALL_DIR accessible through cmake config
- [ ] careful about open mp flag (macro)
- [ ] careful about kamikaze, mpi,
- [ ] careful about march flags
- [ ] Update the scripts for the paraview split again
	> here, generate the create module script at build time and install it under 
	${CMAKE_INSTALL_PREFIX}/bin 
- [x] install everything by default

# Variable naming conventions
- [ ] `TTK_LIBRARY`
- [ ] `TTK_VTK_LIBRARY`
- [ ] `TTK_PARAVIEW_PLUGINS`
- [ ] `TTK_STANDALONE_APPS`
- [x] `TTK_ENABLE_CPU_OPTIMIZATION`
- [x] `TTK_ENABLE_KAMIKAZE`
- [x] `TTK_ENABLE_MPI`
- [x] `TTK_ENABLE_OPENMP`
- [x] `TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE`
- [ ] `TTK_INSTALL_PLUGIN_DIR`
- [ ] change all feature-enabler flags to TTK_ENABLE_XX

# Documentation
- [ ] provide a flag "TTK_BUILD_DOC" that calls doxygen if present, upon install ${CMAKE_INSTALL_PREFIX}/share/doc/ttk/

# ParaView Plugins
- [x] move cmake code to the paraview dir
- [x] fix paraview plugin segfault issues

# Binaries
- [x] move them to an executable location (/usr/local/bin is by default)

# Root CMakeLists.txt
- [-] list all buildadble components (such that users can simply comment things out)
	> there's pros and cons about that feature (let's see that later)
	> let's not do that for now

# third party code
- [x] create a TTK_LIBRARIES variable (just like VTK_LIBRARIES) to make the CMakeLists.txt simpler (see my comment on topo-vol commit)
	> Will: this is done with the bundle alias libraries.

# Travis CI support?
- [x] Build basic TTK libraries
- [ ] Build VTK and ParaView code? This is tough b/c Travis migrated their infrastructure and now we run out of time. We could build the binaries and pull from some server to the CI builder to save this time.
