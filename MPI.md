## Running TTK and Paraview with MPI

Please note that not all TTK filters have MPI support. If a filter does not have MPI support, an error message will be printed at runtime and its behavior will be undefined.

### Compilation

In order to use MPI, Paraview needs to be compiled with MPI support. This can be done by setting the `PARAVIEW_USE_MPI` variable to `ON`. For more information, see [Paraview's building documentation](https://gitlab.kitware.com/paraview/paraview/-/blob/master/Documentation/dev/build.md#linux). TTK also needs to be compiled with MPI support. This can be done by setting the `TTK_ENABLE_MPI` variable to `ON` during compilation. TTK is 
parallelized using both threads and MPI processes simultaneously. Some filters require the use of threads to function using MPI. To measure the computation time of a TTK filter, the variable `TTK_ENABLE_MPI_TIME` needs to be set to `ON`.

### Environment variables
Some TTK filters require a level of thread support not provided by Paraview by default. Obtaining the right level of thread support is MPI implementation dependent. For OpenMPI, it is done by setting the environment variable ` OMPI_MPI_THREAD_LEVEL` to 1. For MPICH, it is done by setting the environment variable `MPIR_CVAR_DEFAULT_THREAD_LEVEL` to 1. 

To specify the number of threads to use during execution, the environment variable `OMP_NUM_THREADS` can be used.

### Execution

To use Paraview in distributed mode using MPI, one has to use either `pvserver` or `pvbatch`. 

Please note that TTK only supports a number of processes equal to $2^n$, with $n\in \mathbb{N}$. There are no constraints on the number of threads.

The use of the periodic grid in distributed mode requires a number of processes equal to $8^n$, with $n\in \mathbb{N}$. It may work for a number equal to $2^n$ if all processes distribute the data similarly in all directions (this distribution depends on Paraview and the dimensions of the data set). 

#### pvserver
The following command allows one to start `pvserver` with 4 MPI processes and 8 threads (per process) with OpenMPI:

    OMPI_MPI_THREAD_LEVEL=1 OMP_NUM_THREADS=8 mpirun -np 4 pvserver

The next step is to run the Paraview client using the command:

    paraview

Now all that is needed is to connect to the server from the client GUI by going into File > Connect. For more information on how to create and connect to a remote Paraview server, see [Section 6.2](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#remote-visualization-in-paraview) and [Section 6.3](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#connect-to-the-remote-server) of Paraview's documentation. For more information on how to use `pvserver` and how it works, see [pvserver's documentation here](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#parallel-processing-in-paraview-and-pvpython).

####pvbatch

`pvbatch` does not require to set up a `pvserver`. However, it does require that the Paraview pipeline to be executed is given in the format of a Python script. To execute the pipeline `pipeline.py` using `pvbatch` using 4 MPI processes and 8 threads (per process) with OpenMPI, the following command can be used:

    OMPI_MPI_THREAD_LEVEL=1 OMP_NUM_THREADS=8 mpirun -n 4 pvbatch pipeline.py

For more information on `pvbatch`, see [pvbatch's documentation here.](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#using-pvbatch)


### Usage within a pipeline

##### For unstructured grids

In order to obtain results identical regardless of the number of processes, the data needs to be prepared by using the filter `TTKIdentifiers` on the data before the distribution. This means opening the data file **in parallel with MPI**, using the `TTKIdentifiers` and then redistributing using the `RedistributeDataSet` filter. The pipeline should look like this:

    Data File > TTKIdentifiers > RedistributeDataSet

After this, the user can start using TTK filters, the rest of the needed preconditioning will be triggered automatically.

##### For structured grids

There is **nothing** to do for structured grids! Everything needed is triggered automatically in the preconditioning steps.

### Output format

- The output of TTK filters executed in distributed mode using MPI is usable as-is in a distributed pipeline.

- The distributed output of a TTK filter is strictly identical, regardless of how many processes are used during execution.
