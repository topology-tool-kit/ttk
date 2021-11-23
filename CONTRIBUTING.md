Thanks for contributing to TTK!

Please find below a few guidelines that we invite you to consider before making a pull request.

# 0. Getting set up with Github
Please find below generic recommendations for setting up your fork of TTK's main repository.
  - Setting up your Github account:
    - Create an account on [Github](https://github.com/).
    - If applicable, we recommand to upgrade (for free) this account to a Github Pro account, through the [Github education program](https://docs.github.com/en/education/explore-the-benefits-of-teaching-and-learning-with-github-education/use-github-in-your-classroom-and-research/apply-for-an-educator-or-researcher-discount).
  - Forking TTK's main repository:
    - Go to [TTK's main source code repository](https://github.com/topology-tool-kit/ttk) and click on the "Fork" button (top right corner).
    - In the remainder, let us designate by `@PUBLIC` the URL of [TTK's main source code repository](https://github.com/topology-tool-kit/ttk) (i.e. `@PUBLIC = https://github.com/topology-tool-kit/ttk`)
    - Similarly, let us designate by `@FORK` the URL of your public fork of TTK's public repository.
  - Creating your private TTK repository:
    - We recommand to use a private repository for the development of unpublished features.
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
    - When developing a new unpublished feature, we recommand to create a new branch on your `@PRIVATE` repository.
    - When this feature is ready to be made public (e.g. after publication of the corresponding research), push the corresponding branch to your `@FORK`. This will enable you to open a pull-request (PR) to the [main TTK repository](https://github.com/topology-tool-kit/ttk).
  - Setting up [ttk-data](https://github.com/topology-tool-kit/ttk-data)
    - The repository [ttk-data](https://github.com/topology-tool-kit/ttk-data) hosts data sets and examples.
    - We recommand that you re-iterate the above procedure (with a public fork and a private repository) for this repository as well.
    - Note that only the features which are covered by examples in [ttk-data](https://github.com/topology-tool-kit/ttk-data) are tested by the continuous integration.

# 1. Authorship
  - Please enter in header files doxygen style information regarding authorship and, if applicable, related publications. See [core/base/topologicalSimplification/TopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/topologicalSimplification/TopologicalSimplification.h) for a base layer example, [core/vtk/ttkTopologicalSimplification/ttkTopologicalSimplification.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkTopologicalSimplification/ttkTopologicalSimplification.h) for a vtk wrapper example and [paraview/TopologicalSimplification/TopologicalSimplification.xml](https://github.com/topology-tool-kit/ttk/blob/dev/paraview/TopologicalSimplification/TopologicalSimplification.xml) for a ParaView plugin example.

# 2. Code formatting
  - To make TTK's source code more homogeneous and readable, we use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) (under Ubuntu, install the package <code>clang-format-11</code>). A style file is already available in TTK's source tree.
**Before creating a new pull request**, please make sure that you clang-formatted your local 
repository by entering the following command at the top of TTK's source tree: <code>
$ clang-format -i -style=file core/\*/\*/\*h core/\*/\*/\*hpp core/\*/\*/\*cpp core/\*/\*/\*inl standalone/\*/\*/\*cpp standalone/\*/\*/\*h standalone/\*/\*cpp </code>

  - To make your life even easier, we recommend that you setup a clang-format pre-commit hook, which will automatically run clang-format on any of your commits to your local repository.
For this, we recommend to use scripts such as [this one](https://github.com/barisione/clang-format-hooks/).

# 3. Examples
- The repository [ttk-data](https://github.com/topology-tool-kit/ttk-data) hosts a list of data sets and example pipelines (primarly stored as ParaView `pvsm` state files). This repository has multiple purposes:
  - Each entry is used to produce a screenshot for the [Gallery page of TTK's website](https://topology-tool-kit.github.io/gallery.html)
  - Each entry serves as a reproducible example with ParaView, as documented on the [Tutorial page of TTK's website](https://topology-tool-kit.github.io/tutorials.html)
  - Each entry is automatically tested by our continuous integration. At the moment, only code features which are used in [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s state files are automatically tested upon pull requests.
  - Each entry is described in detail in [TTK's Example website](https://topology-tool-kit.github.io/examples), for novice users who want to get started with TTK and Python.
- If you develop some new features in TTK (either by creating a new module or extending an existing one), we strongly invite you to produce an entry in the [ttk-data](https://github.com/topology-tool-kit/ttk-data) repository (by pull request). Each entry should be provided as follows:
  - A ParaView `pvsm` state file, stored in the `states` directory:
    - Modify any full path to the input files to make them relative.
  - A Python script (with the same name but with the `py` extension instead of the `pvsm` extension), stored in the `python` directory:
    - This script can be automatically generated by ParaView:
      - Once the state file is opened in ParaView, click on `File`, `Save State...` and make sure to select the entry `Python state file (*.py)` from the extension drop down menu.
    - Edit the automatically generated (and verbose) script:
      - Insert the following line at the top: `#!/usr/bin/env python`
      - After the line `from paraview.simple import *`, remove all the lines *before* and *after* the section entitled ```#setup the data processing pipelines``` (see the automatically generated comments). The main idea is to provide a minimalist and simple Python script which only includes the loading of the input data and the key steps of the data analysis pipeline.
      - Insert new lines at the end of the script to store the outputs of the pipeline with `SaveData()`.
  - A MkDocs file (with the same name as the `pvsm` state file and the `py` Python script) in the `doc` directory. Please create a new entry by copying an already existing one. Each new entry should be organized as follows:
    -
    - Note that the output webpage can be visualized locally by entering the command `mkdocs serve` in the [ttk-data](https://github.com/topology-tool-kit/ttk-data) directory.
  - In the cases (`pvsm` state file, `py` script, MkDocs entry), we invite you to checkout the other examples already included in [ttk-data](https://github.com/topology-tool-kit/ttk-data) for inspiration.

# 4. Code documentation
  -
  - TTK uses [Doxygen](https://www.doxygen.nl/index.html) for the automatic generation of online documentation.
  - Please use Doxygen commands for documenting your header files:
    - Examples for the base layer:
      - [Triangulation.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/triangulation/Triangulation.h)
      - [HelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/helloWorld/HelloWorld.h)
    - Examples for the VTK layer:
      - [ttkHelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkHelloWorld/ttkHelloWorld.h)

# 5. Continuous integration
  - TTK uses some basic continuous integration, which consists in testing for build and run success under Linux, Windows and MacOs upon each commit or pull request. **Your pull request will not be merged if it fails these tests**.


# 6. Submitting code
  - If you plan to submit a **new module**, we invite you to read our [Guidelines for Developing a New TTK Module](https://github.com/topology-tool-kit/ttk/wiki/Guidelines-for-Developing-a-New-TTK-Module). 
  - Prepare your pull-request to the **dev** branch of [TTK](https://github.com/topology-tool-kit/ttk/tree/dev). **Before** submitting it, please make sure that your fork is in sync with the latest version of TTK's source tree (typically by entering a command like <code>git pull ttk-public dev</code>, where <code>ttk-public</code> is the name of your remote pointing to TTK's public source tree). Please make sure that your new code runs fine with TTK's performance mode turned on <code>TTK\_ENABLE\_KAMIKAZE=ON</code> (OFF by default on the **dev** branch).
  - Please submit a pull-request with an example to the **dev** branch of [ttk-data](https://github.com/topology-tool-kit/ttk-data/tree/dev):
    - Provide a ParaView state file (*.pvsm) in the [<code>states/</code>](https://github.com/topology-tool-kit/ttk-data/tree/dev/states) directory which runs your new module. This example will be used to both test and demo your module.
    - Provide new data sets if needed. For new data sets, please update the [README](https://github.com/topology-tool-kit/ttk-data/blob/dev/README) file to include the corresponding references.
    - Make sure to remove all absolute paths from the state file example generated by ParaView (with a text editor). The state files are supposed to be run with ParaView from the root directory of [ttk-data](https://github.com/topology-tool-kit/ttk-data) (see [TTK's tutorial page](https://topology-tool-kit.github.io/tutorials.html) for examples).
  - Please submit a pull-request to [topology-tool-kit.github.io](https://github.com/topology-tool-kit/topology-tool-kit.github.io):
    - For this, generate a screenshot for the above ParaView state file example. Please include on the left of your screenshot a terminal with TTK's output (see examples on [TTK's gallery page](https://topology-tool-kit.github.io/gallery.html)). Please use a dark color theme if possible.
    - Place that screenshot in [<code>img/gallery/</code>](https://github.com/topology-tool-kit/topology-tool-kit.github.io/tree/master/img/gallery)
  - Prepare, if possible, a video tutorial, similar to those available on [TTK's tutorial page](https://topology-tool-kit.github.io/tutorials.html). In this video, you should:
    - Open a terminal and load your state file demo to explain what it does and what it shows.
    - Re-open ParaView and re-create your state file pipeline **from scratch** (to show people how to put things together)

From our experience, video tutorials are **essential** to attract users towards new features. So if you want people to use your new module, you want to prepare a video tutorial. There are several excellent software packages for video editing (for instance [kdenlive](https://kdenlive.org/en/))
