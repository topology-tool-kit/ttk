Thanks for contributing to TTK!

Please find below a few guidelines that we invite you to consider before making a pull request.

# 0. Getting set up with Github
Please find below generic recommendations for setting up your fork of TTK's main repository.
  - Setting up your Github account:
    - Create an account on [Github](https://github.com/).
    - If applicable, we recommend to upgrade (for free) this account to a Github Pro account, through the [Github education program](https://docs.github.com/en/education/explore-the-benefits-of-teaching-and-learning-with-github-education/use-github-in-your-classroom-and-research/apply-for-an-educator-or-researcher-discount).
  - Forking TTK's main repository:
    - Go to [TTK's main source code repository](https://github.com/topology-tool-kit/ttk) and click on the "Fork" button (top right corner).
    - In the remainder, let us designate by `@PUBLIC` the URL of [TTK's main source code repository](https://github.com/topology-tool-kit/ttk) (i.e. `@PUBLIC = https://github.com/topology-tool-kit/ttk`)
    - Similarly, let us designate by `@FORK` the URL of your public fork of TTK's public repository.
  - Creating your private TTK repository:
    - We recommend to use a private repository for the development of unpublished features.
    - For this, create and setup a *private* repository (e.g. `ttk-yourusername`). Let us call it `@PRIVATE`.
    - Clone your `@PRIVATE` repository locally and enter the following commands:
    ```
    $ git clone @PRIVATE ttk-yourusername
    $ cd ttk-yourusername
    $ git remote add ttk-public @PUBLIC
    $ git remote add ttk-fork @FORK
    ```
  - Daily usage:
    - At this point, you can regularly keep your private repository up-to-date by entering the following command in it:
    ```
    $ git pull ttk-public dev
    ```
    - When developing a new unpublished feature, we recommend to create a new branch on your `@PRIVATE` repository.
    - When this feature is ready to be made public (e.g. after publication of the corresponding research), push the corresponding branch to your `@FORK`. This will enable you to open a pull-request (PR) to the [main TTK repository](https://github.com/topology-tool-kit/ttk).
  - Setting up [ttk-data](https://github.com/topology-tool-kit/ttk-data)
    - The repository [ttk-data](https://github.com/topology-tool-kit/ttk-data) hosts data sets and examples.
    - We recommend that you re-iterate the above procedure (with a public fork and a private repository) for this repository as well.
    - Note that only the features which are covered by examples in [ttk-data](https://github.com/topology-tool-kit/ttk-data) are tested by the continuous integration.

# 1. Authorship
  - Please enter in header files doxygen style information regarding authorship and, if applicable, related publications. See [core/base/topologicalSimplification/TopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/topologicalSimplification/TopologicalSimplification.h) for a base layer example, [core/vtk/ttkTopologicalSimplification/ttkTopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkTopologicalSimplification/ttkTopologicalSimplification.h) for a vtk wrapper example and [paraview/TopologicalSimplification/TopologicalSimplification.xml](https://github.com/topology-tool-kit/ttk/blob/dev/paraview/TopologicalSimplification/TopologicalSimplification.xml) for a ParaView plugin example.

# 2. Code formatting
  - To make TTK's source code more homogeneous and readable, we use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) (under Ubuntu, install the package <code>clang-format-11</code>). A style file is already available in TTK's source tree.
**Before creating a new pull request**, please make sure that you clang-formatted your local
repository by entering the following command at the top of TTK's source tree:

```sh
$ git ls-files | grep -E "\.cpp$|\.cxx$|\.h$|\.hpp$|\.inl$" | xargs -P$(nproc) -n1 clang-format -i
```

  - `clang-format` can also be used on individual files (`$
    clang-format -i source.cpp`) or on Git diffs with
    `git-clang-format`.

  - To make your life even easier, we recommend that you setup a clang-format pre-commit hook, which will automatically run clang-format on any of your commits to your local repository.
For this, we recommend to use scripts such as [this one](https://github.com/barisione/clang-format-hooks/).

# 3. Examples
- The repository [ttk-data](https://github.com/topology-tool-kit/ttk-data) hosts a list of data sets and example pipelines.
- If you develop some new feature in TTK (either by creating a new module or by extending an existing one), we strongly invite you to produce an entry in the [ttk-data](https://github.com/topology-tool-kit/ttk-data) repository (by pull request).
- Note that, at the moment, only code features which are used in [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s state files are automatically tested by continuous integration (upon code pull requests). Also, note that only such code features (i.e. automatically tested by continuous integration) are considered for integration in Kitware's ParaView official distribution.
- Please checkout [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s [CONTRIBUTOR Guide](https://github.com/topology-tool-kit/ttk-data/blob/dev/CONTRIBUTING.md) for detailed instructions.

# 4. Code documentation
  - TTK uses [Doxygen](https://www.doxygen.nl/index.html) for the automatic generation of online documentation.
  - Each base layer header file should be organized as follows:
    - A overview, including in order:
      - the line `\ingroup base`
      - the name of the class, with the command `\class ttk::`
      - the author(s) of the class, with the command `\author`
      - the date of creation of the class, with the command `\date`
      - a one-liner description of the class, with the command `\brief`
      - a long description of the purpose of the class
      - related publications (if applicable)
      - a list of related classes, with the command `\sa`:
        - other base layer modules commonly used in conjunction with the present filter
      - a pointer to the corresponding VTK layer class (if applicable), with the command `\sa`
    - Each function should be documented as follows:
      - a quick description, with preconditions `\pre`, notes `\note` or warnings `\warning` if needed.
      - a description of each attribute, with the command `\param`
      - a description of the return code, with the command `\return`
      - a list of related functions, with the command `\sa`
  - Each VTK layer header file should be organized as follows:
    - A overview, including in order:
      - the line `\ingroup vtk`
      - the name of the class, with the command `\class`
      - the author(s) of the class, with the command `\author`
      - the date of creation of the class, with the command `\date`
      - a one-liner description of the class, with the command `\brief`
      - a long description of the purpose of the class
      - if the class is a filter (this is the case in general):
        - a description of each input, with the command `\param`, describing:
          - its class (e.g. `vtkDataSet`, `vtkPointSet`, etc.)
          - the required data arrays if any (type, size, meaning).
        - a description of each output, with the command `\param`, describing:
          - the output class (e.g. `vtkDataSet`, `vtkPointSet`, etc.)
          - a description of each new data array:
            - type, size, meaning
            - if a new data are encodes with an integer code a specific classification (for instance, a region type), please clarify the mapping class/integer (e.g. the correspondence integer to region type).
        - recall that the class can be used as any other VTK filter
      - related publications (if applicable)
      - a list of related classes, with the command `\sa`:
        - other VTK filters commonly used in conjunction with the present filter
      - a pointer to the corresponding base code class (if applicable), with the command `\sa`
      - the list of **all** the entries from [TTK's Examples website](https://topology-tool-kit.github.io/examples/) including the filter.
  - Each ParaView XML file should be organized as follows:
    - an overview `<Documentation>` section, within the main `<SourceProxy>` section:
      - this should be a copy of the overview description of the corresponding VTK class (see above).
    - for each `<InputProperty>` section:
      - include a `<Documentation>` section describing the corresponding input data, by specifying:
        - its class (e.g. `vtkDataSet`, `vtkPointSet`, etc.)
        - the required data arrays if any (type, size, meaning).
    - for each other Property section (`<IntProperty>`, `<StringProperty>`, etc.), specify:
      - the meaning of the input option and its possible values.
  - Please check out the following examples for inspiration:
    - Examples for the base layer:
      - [Triangulation.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/triangulation/Triangulation.h)
      - [TopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/topologicalSimplification/TopologicalSimplification.h)
      - [HelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/helloWorld/HelloWorld.h)
    - Examples for the VTK layer:
      - [ttkTopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkTopologicalSimplification/ttkTopologicalSimplification.h)
      - [ttkHelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkHelloWorld/ttkHelloWorld.h)
    - Examples for the ParaView XML layer:
      - [TopologicalSimplification](https://github.com/topology-tool-kit/ttk/blob/dev/paraview/xmls/TopologicalSimplification.xml)


# 5. Continuous integration

- TTK uses some basic continuous integration, which consists in
  testing for build and run success under Linux, Windows and MacOs
  upon each commit or pull request. **Your pull request will not be
  merged if it fails these tests**.

- TTK developers can either try to launch themselves the CI workflows
  on their public fork (see [the CI
  documentation](./.github/workflows/README.md)) or replicate locally
  the checks performed by the CI (way faster).

- The CI makes heavy use of Clang tools, such as `clang-check` or
  `clang-tidy`. The first step is to use CMake to generate a
  Compilation Database, a JSON file that contains explicit compilation
  commands for each source file in the repository. Go to your build
  directory, then

  ```sh
  $ cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE .
  ```

  This will generate a `compile_commands.json` file in your build
  directory.

- The `clang-check` tool parses source files according to the
  corresponding compilation database entry and displays eventual
  compilation errors or warnings. This tool can be used to complement
  GCC (it is quicker to run and the error descriptions are usually
  better phrased). From TTK's root directory, use `$ clang-check -p
  build source.cpp` where `build` is the build directory (which
  contains the compilation database). To pass over all TTK source
  files, use:

  ```sh
  $ git ls-files | grep core | grep -E "\.cpp$|\.cxx$" | xargs -P$(nproc) -n1 clang-check -p build
  ```

- Among the warnings that are tested with `clang-tidy` in the CI are
  the `unused-parameter` warnings. To silence them, TTK offers several solutions:

  * the `ttkNotUsed` macro used on the parameter declaration
    highlights the non-usage of the parameter in the function body;
  * the `TTK_FORCE_USE` macro in the function body is treated by the
    compiler as a parameter use.

  See the definition of these variables in
  [BaseClass.h](./core/base/common/BaseClass.h) for use-cases.

- Finally, the `clang-tidy` tool is used in the CI to perform C++ code
  style checks and some basic static analysis (using Clang's Static
  Analyzer). The checks are listed in the [.clang-tidy](./.clang-tidy)
  file. A description of each of these checks is available on
  [clang-tidy's
  website](https://clang.llvm.org/extra/clang-tidy/checks/list.html).
  Similarly to `clang-check`, `clang-tidy` can be used on individual
  source files with `$ clang-tidy -p build source.cpp`, or on the
  whole TTK codebase:

  ```sh
  $ git ls-files | grep core | grep -E "\.cpp$|\.cxx$" | xargs -P$(nproc) -n1 clang-tidy -p build
  ```

  Note that a `--fix` flag is available to automatically correct some
  warnings.

- Static analysis is performed using the
  `--checks="-*,clang-analyzer-*"` flags that disables every
  `clang-tidy` warning except the Clang Static Analyzer ones:

  ```sh
  $ git ls-files | grep core | grep -E "\.cpp$|\.cxx$" \
    | xargs -P$(nproc) -n1 clang-tidy --checks="-*,clang-analyzer-*" -p build
  ```

- Code editors, thanks to the [Language Server
  Protocol](https://microsoft.github.io/language-server-protocol/),
  can take advantage of `clang-check` and `clang-tidy` and display
  their warnings in their interface.

# 6. Submitting code
  - If you plan to submit a **new module**, we invite you to read our [Guidelines for Developing a New TTK Module](https://github.com/topology-tool-kit/ttk/wiki/Guidelines-for-Developing-a-New-TTK-Module).
  - Prepare your pull-request to the **dev** branch of [TTK](https://github.com/topology-tool-kit/ttk/tree/dev). **Before** submitting it, please make sure that your fork is in sync with the latest version of TTK's source tree (typically by entering a command like <code>git pull ttk-public dev</code>, where <code>ttk-public</code> is the name of your remote pointing to TTK's public source tree). Please make sure that your new code runs fine with TTK's performance mode turned on <code>TTK\_ENABLE\_KAMIKAZE=ON</code> (OFF by default on the **dev** branch).
  - Please submit a pull-request with an example to [ttk-data](https://github.com/topology-tool-kit/ttk-data/tree/dev). See [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s [CONTRIBUTOR Guide](https://github.com/topology-tool-kit/ttk-data/blob/dev/CONTRIBUTING.md) for detailed instructions.
  - Prepare, if possible, a video tutorial, similar to those available on [TTK's tutorial page](https://topology-tool-kit.github.io/tutorials.html). In this video, you should:
    - Open a terminal and load your state file demo to explain what it does and what it shows.
    - Re-open ParaView and re-create your state file pipeline **from scratch** (to show people how to put things together)

From our experience, video tutorials are **essential** to attract users towards new features. So if you want people to use your new module, you want to prepare a video tutorial. There are several excellent software packages for video editing (for instance [kdenlive](https://kdenlive.org/en/))

# 7. Tools for TTK developers

A list of tips & tools to help TTK developers has been gathered into
[DEVELOPER_EXPERIENCE.md](./doc/DEVELOPER_EXPERIENCE.md).
