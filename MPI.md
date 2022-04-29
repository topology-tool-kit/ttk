## Running TTK and Paraview with MPI

### Compilation

In order to use MPI, Paraview needs to be compiled with MPI support. This can be done by setting the `PARAVIEW_USE_MPI` variable to `ON`. For more information, see [Paraview's building documentation](https://gitlab.kitware.com/paraview/paraview/-/blob/master/Documentation/dev/build.md#linux). TTK also needs to be compiled with MPI support. This can be done by setting the `TTK_ENABLE_MPI` variable to `ON` during compilation. TTK is 
parallelized using both threads and MPI processes simultaneously. Some filters require the use of threads to function using MPI.

### Environment variables
Some TTK filters require a level of thread support not provided by Paraview by default. Obtaining the right level of thread support is MPI implementation dependent. For OpenMPI, it is done by setting the environment variable ` OMPI_MPI_THREAD_LEVEL` to 3. For MPICH, it is done by setting the environment variable `MPIR_CVAR_DEFAULT_THREAD_LEVEL` to 3. 

To specify the number of threads to use during execution, the environment variable `OMP_NUM_THREADS` can be used.

### Execution

To use Paraview in distributed mode using MPI, one has to use either `pvserver` or `pvbatch`. 

#### pvserver
The following command allows one to start `pvserver` with 4 MPI processes and 8 threads (per process) with OpenMPI:

    OMPI_MPI_THREAD_LEVEL=3 OMP_NUM_THREADS=8 mpirun -np 4 pvserver

The next step is to run the Paraview client using the command:

    paraview

Now all that is needed is to connect to the server from the client GUI by going into File > Connect. For more information on how to create and connect to a remote Paraview server, see [Section 6.2](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#remote-visualization-in-paraview) and [Section 6.3](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#connect-to-the-remote-server) of Paraview's documentation. For more information on how to use `pvserver` and how it works, see [pvserver's documentation here](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#parallel-processing-in-paraview-and-pvpython).

####pvbatch

`pvbatch` does not require to set up a `pvserver`. However, it does require that the Paraview pipeline to be executed is given in the format of a Python script. To execute the pipeline `pipeline.py` using `pvbatch` using 4 MPI processes and 8 threads (per process) with OpenMPI, the following command can be used:

    OMPI_MPI_THREAD_LEVEL=3 OMP_NUM_THREADS=8 mpirun -n 4 pvbatch pipeline.py

For more information on `pvbatch`, see [pvbatch's documentation here](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#using-pvbatch)


### Usage within a pipeline

#### Before distribution of the data

In order to obtain results identical to results obtained if the pipeline is executed in sequential, the data needs to be prepared by using the filter `GenerateGlobalIds` on the data in sequential. The output of this filter should then be used for the distributed computation.

#### After distribution of the data

 - TTK filters require ghost cells and ghost points to be created to function correctly. Ghost cells are cells inside a data set that are copies of the interfacing cells of an adjacent data set. TTK filters require two layers of ghost cells. Ghost cells and points can be created using the [`GhostCellsGenerator` filter](https://kitware.github.io/paraview-docs/latest/python/paraview.simple.GhostCellsGenerator.html).

 - The `GhostCellPreconditioning` filter can be used after the  `GhostCellsGenerator` filter to calculate which rank is the "primary" owner of each vertex.

- The `ArrayPreconditioning` filter requires both of the previous filters to function correctly in distributed mode.

### Output format

- The output of TTK filters executed in distributed mode using MPI is usable as-is in a distributed pipeline. There is no need to add ghost cells or global ids to the output of a TTK filter.

- The distributed output of a TTK filter is strictly identical, regardless of how many processes are used during execution.
