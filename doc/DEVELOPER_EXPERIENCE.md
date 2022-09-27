TTK Developer Experience
========================

Author: Pierre Guillou
Date: December 2022

This document contains a list of tips and pieces of advice that can
help TTK developers to increase their effectiveness.

## Software of interest

Here is a list of software applications of interest:

* [CMake][1] is the build tool used by TTK. By default, it generates
  Makefile recipes.
* [Ninja][2] is a faster and more modern alternative to Makefiles. It
  is parallel by default.
* [clang-format][3] is the tool used to format TTK's C++ code.
* [clang-check][4] allows to quickly check for compiler warnings or
  errors with Clang messages (often more helpful than GCC ones).
* [clang-tidy][5] performs extended checks on C++ code.
* [clangd][6] is a [LSP server][13] that integrates `clang-format` and `clang-check`.
* [ccache][7] accelerates the build process by storing object files
  into a separate cache.
* [perf][10] is a Linux tool that collects performance statistics at
  run-time.
* [hotspot][11] helps at visualizing `perf` traces.
* [heaptrack][12] collects and displays the run-time memory footprint.

## Environment variables

These environment variables are used by [CMake][1] to customize the
build process:

* `CMAKE_GENERATOR=Ninja` replaces GNU Make by [Ninja][2] as build
  tool.
* `CMAKE_COLOR_DIAGNOSTICS=ON` enables coloration on compiler warnings
  with Ninja (coloration is enabled by default with GNU Make).
* `CMAKE_EXPORT_COMPILE_COMMANDS=ON` generates a [compilation
  database][8] in the build directory. This compilation database can
  be fed to useful tools such as [`clang-check`][4], [`clang-tidy`][5]
  or [`clangd`][6].
* `CMAKE_CC_COMPILER_LAUNCHER=ccache` and
  `CMAKE_CXX_COMPILER_LAUNCHER=ccache` make CMake use [`ccache`][7] to
  store generated binary objects in a cache to make faster builds
  after clean.

## TTK CMake keys

These TTK-specific CMake keys can be tweaked to customize the TTK
build.

* `CMAKE_CXX_FLAGS` set `-ggdb -fno-omit-frame-pointer` adds debug
  data into TTK libraries. A debugger (such as GDB) needs it to
  provide accurate information.
* `TTK_TIME_TARGETS` prints the time taken by the compilation of every
  translation unit (`.cpp` file). Quite useful for learning where the
  build time is spent.
* `TTK_REDUCE_TEMPLATE_INSTANTIATIONS` uses a reduced list of template
  instantiations for the VTK layer. May not work with all `ttk-data`
  state files.
* `TTK_SCRIPTS_PATH` is only used by the `DimensionReduction` module
  by now. It points to the folder that should contain the
  `dimensionReduction.py` Python script file (at run-time). By
  default, it points to the install prefix but it can be modified to
  point to the relevant source directory (to skip the installation process).

## Format C++ code

Each time a PR is opened on GitHub, the CI checks that the provided
code is correctly formatted using [`clang-format`][3]. The formatting
step is often forgotten in the PR, so here is a one-liner that uses
Git to format everything that have changed since the last commit in
the reference repository (here the `origin` remote).

`git-clang-format origin/dev`

## Better, more understandable warnings

[clang-check][4] is a CLang tool that quickly checks the syntax of
C/C++ code without generating assembly to provide warnings that often
are more understandable that GCC ones. To use it, generate a
compilation database (`CMAKE_EXPORT_COMPILE_COMMANDS=ON`) in the build
directory, then call `clang-check` on the failing C++ source file (`.cpp`):

```bash
clang-check -p build core/vtk/ttkPersistenceDiagram.cpp
```

Here, the `-p build` passes the path to the build directory to
`clang-check`.

[clang-tidy][5] is quite similar but offers more extensive checks of
C/C++ code but takes more time to complete. It is use the same way as
`clang-check`. There's a list of enabled and disabled warnings in the
[.clang-tidy](../.clang-tidy) file.

## Launch CI pipelines

CI workflows can be launched in GitHub forks (here the `public`
remote) of the TTK repository by pushing tags.

For launching the `check_code` workflow, use

```
git tag -f check && git push -f public check
```
For launching the `test_build` workflow, use

```
git tag -f test_build && git push -f public test_build
```
For launching the `packaging` workflow, use

```
git tag -f package && git push -f public package
```

## Install ParaView in a separate prefix

To avoid cluttering the `/usr/local` prefix, I advise to install
ParaView in a prefix inside the user's `$HOME`. To do so, first set
the `CMAKE_INSTALL_PREFIX` to the desired prefix (for instance
`$HOME/install`) during ParaView configuration. Once ParaView is built
and installed inside this prefix, setting these environment variables
makes ParaView available to the shell.

```bash
PV_PREFIX=$HOME/install/
export PATH=$PV_PREFIX/bin:$PATH
export CMAKE_PREFIX_PATH=$PV_PREFIX/lib/cmake:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PV_PREFIX/lib/python<ver>/site-packages:$PYTHONPATH
```

## Use TTK from build directory (no install)

Similarly, TTK can also be installed in a separate prefix OR we can
skip the installation process and use directly the binaries located
inside the build directory (here `$HOME/ttk/build`). These environment
variables helps ParaView to find the TTK plugin:

```bash
TTK_BUILD_DIR=$HOME/ttk/build
export PATH=TTK_BUILD_DIR/bin:$PATH
export LD_LIBRARY_PATH=$TTK_BUILD_DIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$TTK_BUILD_DIR/lib/python<ver>/site-packages:$PYTHONPATH
export PV_PLUGIN_PATH=$TTK_BUILD_DIR/lib/TopologyToolKit:$PV_PLUGIN_PATH
```

The `TTK_SCRIPTS_PATH` CMake key should also point to the source
directory `core/base/dimensionReduction/` to find the
`dimensionReduction.py` Python script at run-time.

## Debug with GDB

[GDB][14] is the default debugger of the GNU Project. It it mostly
used to stop the execution of a program at any point to examine all
the local variables. For TTK, it is possible to use GDB at the module
level by including debugging information and reducing the compiler
optimization level. Edit the local `CMakeLists.txt` file inside the
TTK <module> directory and add this line:

```cmake
target_compile_options(<module> PUBLIC -ggdb -O0)
```

GDB can now track individual variables in the program.

GDB supports breakpoints, i.e. stopping the execution of the program
on a specific line of code. However, it is quite hard to use this
feature with TTK when used through ParaView, as TTK is dynamically
loaded during the ParaView execution. Besides, GDB breakpoints often
slow down the execution. Fortunately, there is a way to insert a
breakpoint inside the C++ source code: just use the `kill` C function
with the `SIGINT` signal, as displayed below.

```c
// for SIGINT
#include <csignal>

int foo() {
  kill(getpid(), SIGINT);
}
```

The debugger will stop the execution at the `kill` line and give back
control to the user, who can now examine the variables and continue
the execution, either line by line or untile the termination of the
program.

## Debug with the Address Sanitizer

The [Address Sanitizer][9] tracks memory accesses at run-time and
helps solving segfaults. Inside TTK, it can be enabled for individual
modules. In the `CMakeLists.txt` file inside the TTK <module>
directory, add these lines:

```cmake
target_compile_options(<module> PUBLIC -fsanitize=address)
target_link_options(<module> PUBLIC -fsanitize=address)
```

The `LD_PRELOAD` environment variable should be set to preload the
ASan libraries at run-time.

The Address Sanitizer really shines in conjunction with GDB. With
those lines in `$HOME/.gdbinit`:

```
set environment ASAN_OPTIONS=abort_on_error=1:detect_leaks=0
set environment LD_PRELOAD=/usr/lib/libasan.so
```

GDB will be able to stop at the first ASan error. The user then can
examine the backtrace. Note that ParaView takes a very long time to
load with ASan pre-loaded. It is often more convenient to debug a
Python script launched by the Python interpreter than a ParaView state
file.

Launch your program inside GDB with `gdb -ex run --args <your
command>`.

## Avoid manual allocations

`new`, `delete` and their C counterparts `malloc` and `free` are not
needed anymore in modern C++. They can easily be replaced by
`std::vector`, `std::shared_ptr` or `std::unique_ptr` which prevent
memory leaks. Sometimes the heap object can even be placed on the
stack. TTK is currently mostly free of manual memory management,
please keep it this way.

## Visualize performance & memory footprint

To better know where a TTK module spends computation time, the
[perf][10] tool can be used to collect traces from performance
counters on an execution of any module.

For instance, in `ttk-data`, the command

```bash
perf record --call-graph dwarf paraview --state=states/dragon.pvsm
```

records in the `perf.data` file the traces of the execution of the
`dragon.pvsm` state file (also working with Python scripts and
standalones).

If traces seem to be missing data, don't forget to build TTK with
debugging info (`-ggdb -fno-omit-frame-pointer`).

To visualize the traces, call [hotspot][11] directly inside the
directory that contains the `perf.data`. Hotspot provides interesting
visualizations of how the CPU cycles are spent with the Flamegraph.

To know the memory footprint, use [heaptrack][12]. This time, the
traces generator and the visualization application are the same.

```bash
heaptrack paraview --state=states/dragon.pvsm
```

will record the memory consumption of the state file, then open the
visualization application.

## Work with one local TTK Git repository

`git remote` allow a local Git repository to point to several
remotes. It can be useful for fetching TTK latest commits from the
reference repository (https://github.com/topology-tool-kit/ttk) and
pushing them into a GitHub fork.

To do so, use

```bash
git remote add <remote-name> git@github:<user>/<repo>
```

to add another remote to the local repository.

Then, `git fetch --all` will synchronize & fetch the latest commits
from GitHub into the local repository. `git pull <remote> <branch>`
will synchronize a local branch with its remote.

[1]: https://cmake.org/
[2]: https://ninja-build.org/
[3]: https://clang.llvm.org/docs/ClangFormat.html
[4]: https://clang.llvm.org/docs/ClangCheck.html
[5]: https://clang.llvm.org/extra/clang-tidy/
[6]: https://clangd.llvm.org/
[7]: https://ccache.dev/
[8]: https://clang.llvm.org/docs/JSONCompilationDatabase.html
[9]: https://github.com/google/sanitizers/wiki/AddressSanitizer
[10]: https://perf.wiki.kernel.org/index.php/Main_Page
[11]: https://github.com/KDAB/hotspot
[12]: https://github.com/KDE/heaptrack
[13]: https://microsoft.github.io/language-server-protocol
[14]: https://www.sourceware.org/gdb/
