TTK's Continuous Integration
============================

TTK uses Continuous integration (CI) each time a contributor submits
new code to TTK, through GitHub's Pull Requests (PR) mechanism.

The CI performs a variety of checks on the last commit of the PR
branch that try to guarantee some static and run-time properties of
the proposed code (no new warnings, no segfaults, etc.).

Currently, TTK's CI relies on the GitHub Actions platform. GitHub
Actions parses and executes Workflow files contained in the
[.github/workflows](../.github/workflows) directory.

Additionally, the [.github/actions](../.github/actions) folder defines
Composite Actions that are reused inside Workflow files.

Some principles I wanted to follow:
* The least surprise: Don't launch workflows in TTK forks
* Efficiency: Avoid doing the same thing twice
* User-friendly: Users should be able to launch the CI by themselves
* Fastness: CI should complete quickly
* Extensiveness: CI should test everything testable
* Correctness: Breaking changes should be highlighted

[GitHub's documentation](https://docs.github.com/en/actions) is pretty
extensive about these Actions.

Workflows
---------

### YAML

GitHub Actions workflows are written in YAML, a description language
similar to JSON, TOML or XML. One pitfall of writing YAML is
indentation, which is very strict. There are some YAML formatting
tools such as [Prettier](https://prettier.io/) but also a [YAML LSP
server](https://github.com/redhat-developer/yaml-language-server) that
can detect these formatting issues.

There are 3 main GitHub Actions workflows for TTK.

### [`check_ci`](../.github/workflows/check.yml)

The main purpose of this workflow is to ensure that the code quality
stays relatively good. To do so, it leverages [clang
tools](https://clang.llvm.org/docs/ClangTools.html) (and a small set
of custom checks). The workflow is divided into 3 main jobs that use
job matrices.

It can be triggered by pushing the `check` tag.

#### The `check-formatting` job: code formatting

This job uses `clang-format` to check that TTK's C++ code is correctly
formatted. It returns an error if not (it could commit changes but I
kind of dislike auto commits). Other checks ensure consistency in the
codebase: Unix line endings, non-empty files, no VTK includes in the
base layer, no `using namespace` directive in header files (can have
side-effets).

#### The `lint-code` job: linting, static Analysis & Doxygen check

This job uses a matrix to perform these three actions in parallel:

* it uses `clang-tidy` to check the codebase against a list of rules
  defines in the [`.clang-tidy` file](../.clang-tidy);
* it uses `clang-tidy` to perform a very basic static analysis on the
  codebase;
* it checks the Doxygen documentation against the method signatures.

The static analysis step is the longest of the whole workflow (~25
minutes). We might want to execute it only on the modified code (but
finding all source files that include a header file might be hard) or
get rid of it if it takes too long.

#### The `check-warnings` jobs: avoid compiler warnings

This job uses `clang-check` to ensure that no warning is introduced in
TTK (using Clang as a compiler, GCC may have different warnings). By
extension, it also verifies if TTK can be built without errors in a
variety of configurations given by the input job matrix. Currently, 5
binary configurations are mixed in the matrix, amounting to 32
effective jobs running in parallel (the limit should be 128). I think
there's a limit of 16 jobs that can be run concurrently in a workflow,
so about half the jobs are postponed till the first half finishes
(this takes around 5 minutes per job).

### [`packaging`](../.github/workflows/package.yml)

This workflow aims at generating binary packages for TTK that can be
user-installed. To do so, corresponding TTK-ParaView packages need to
be built. This is done with the...

#### [ParaView `packaging`](../paraview/patch/package.yml)

This workflow builds pre-patched Paraview packages with the TTK
branding. These packages should be user-installable.

On macOS and Windows, third-party package managers are used to provide
build- and run-time dependencies for both ParaView and TTK: Homebrew
and Anaconda. It should be easy to use these package managers to
install the required dependencies on users' machines. Since a separate
Python installation is available on the GitHub VMs, these separate
installations should be removed before the CMake configuration step to
correctly target brew's or conda's Python.

For building TTK packages, the workflow is divided in two: first,
platform-dependent jobs build the packages, then new jobs are fired to
test these packages outside the build environment.

On macOS, we have not (yet) found a easy-to-use package
format. Instead, we rely on binary archives (`.tar.gz` of the
installation prefix). As a consequence, there's no clear
uninstallation process.

Since TTK is now embedded inside the official ParaView binaries, there
is less rationale to use these packaging workflows.

### [`test_build`](../.github/workflows/test.yml)

This workflow is a simpler variation of the `packaging` workflow. It
is run on every pull request. It checks if TTK builds correctly on the
three main platforms then runs
[ttk-data](https://github.com/topology-tool-kit/ttk-data) state files
(on Ubuntu) and Python scripts to highlight changes. ParaView binaries
still have to be pre-built for this job. This is done in...

#### [ParaView `build_headless_packages`](../paraview/patch/headless.yml)

This is a simpler variation of the ParaView `packaging` workflow. The
main difference is that the Qt GUI is not built. We use
[OSMesa](https://docs.mesa3d.org/osmesa.html) on Ubuntu to do
Offscreen Rendering.

The Offscreen Rendering feature is needed to render the ttk-data state
files into screenshots that are compared to a screenshot database
stored on the ttk-data repository. This system kind of also works on
macOS, but the rendering is way slower than on Ubuntu, so it is
disabled. On Windows, we have not found a way yet to do Offscreen
Rendering and the renderer segfaults on the headless GitHub servers.

All platform execute the ttk-data Python scripts. These scripts
generate output datasets in the VTK formats. The generated files are
hashes (SHA1) then the hashes are compared to a reference
(per-platform) hashes database, also stored in ttk-data.

Since the workflow builds TTK regularly, `ccache` (`sccache` on
Windows) is used to store object files between runs. GitHub provides a
cache system, however it is not shared between branches or PRs. We
rely instead on GitHub (pre-)releases to hold these binary assets
between PRs. The assets are updated each time the workflow is run on
the `dev` branch of the reference repository (after every PR merge).

Debugging
---------

### Issue following VM update

The GitHub virtual machines are often updated, which can sometimes
breaks our workflows. Usually, there are announcements about major
updates on the dedicated [runner image
webpage](https://github.com/actions/runner-images) on GitHub.

One of the most common issues is a software update in a library shared
by ParaView and TTK (e.g. Python). This mostly happens on macOS and
Windows, Ubuntu tend to pin the versions of their packaged software
per Ubuntu release.

2 workarounds exist:
- either pin the packaged version of the offending software to the one
  used to build the TTK-ParaView binaries (`python@3.10` in `brew
  install`, `python=3.10` in `conda install`)
- or make Julien rebuild the reference TTK-ParaView binaries.

More complex issues are hard to track down since GitHub Actions VMs
are not downloadable. Two methods are possible: first, trying to
reproduce the virtual environment and secondly, trying with a custom
workflow on GitHub Actions. The custom workflow can be extracted from
an existing one and be reduced to a very simple and fast run by
disabling dependencies or removing TTK code that takes a long time to
compute. Once the problem is solved on the custom workflow, the
solution can be applied on existing ones.

### Hashes/screenshots discrepancy

After several updates, the hashes & screenshots generated by the CI
may be out of sync with the reference stored in ttk-data.

For the screenshots, the CI generates an archive (per Ubuntu) of
diverging screenshots. These screenshots can then be pasted into
the local ttk-data repository and PR'd into the GitHub reference.

In the `test_build` workflow, a special step filters the ParaView
state files to be executed and the screenshots to be considered. This
can be useful if one screenshot is not deterministic across CI runs.

For the hashes, the CI outputs the content of the generated hashes
file before comparing it to the reference. The reference can then be
overwritten by the corresponding generated hashes.

To this day, several Python scripts generate non-deterministic
outputs. The offending TTK modules are `DimensionReduction` (which
calls `scikit-learn` internally) and `PersistenceDiagramClustering`
(which uses a progressive algorithm that iterates under a time
constraint). TTK modules should ideally be deterministic.
