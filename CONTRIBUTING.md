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
- The repository [ttk-data](https://github.com/topology-tool-kit/ttk-data) hosts a list of data sets and example pipelines.
- If you develop some new feature in TTK (either by creating a new module or by extending an existing one), we strongly invite you to produce an entry in the [ttk-data](https://github.com/topology-tool-kit/ttk-data) repository (by pull request).
- Note that, at the moment, only code features which are used in [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s state files are automatically tested by continuous integration (upon code pull requests). Also, note that only such code features (i.e. automatically tested by continuous integration) are considered for integration in Kitware's ParaView official distribution.
- Please checkout [ttk-data](https://github.com/topology-tool-kit/ttk-data)'s [CONTRIBUTOR Guide](https://github.com/topology-tool-kit/ttk-data/blob/dev/CONTRIBUTING.md) for detailed instructions.

# 4. Code documentation
  - TTK uses [Doxygen](https://www.doxygen.nl/index.html) for the automatic generation of online documentation.
  - Please use Doxygen commands for documenting your header files:
    - Examples for the base layer:
      - [Triangulation.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/triangulation/Triangulation.h)
      - [HelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/base/helloWorld/HelloWorld.h)
    - Examples for the VTK layer:
      - [ttkHelloWorld.h](https://github.com/topology-tool-kit/ttk/blob/dev/core/vtk/ttkHelloWorld/ttkHelloWorld.h)
  - Each VTK layer header file should be organized as follows
    - A overview, including in order:
      - the line `\ingroup vtk`
  - Each ParaView XML file should embed the same documentation description as provided in the corresponding VTK layer header file.

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
