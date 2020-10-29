# Paraview + TTK Docker Image

This docker image contains an installation of the [Topology Tool Kit (TTK)](http://topology-tool-kit.github.io) and [ParaView](http://www.paraview.org) server built from source. Specifically:

- ParaView server with offscreen rendering using either [OSMesa](http://www.mesa3d.org/osmesa.html) or [OSPRay](http://www.ospray.org).
- TTK for ParaView plugins are installed.

It is designed to be used in conjunction with a local ParaView GUI, but can also be used to execute ParaView python scripts.

## Simple usage

To run Kitware's binary distribution of ParaView with TTK's docker, simply run:

``` ./runParaViewTTKDocker.sh <Path to ParaView binary (supported version)> [<Standard ParaView arguments (state files, data, etc.)>]```

To run a python script which uses TTK, simply run:

``` ./runTTKPythonDocker.sh [<Standard pvpython arguments: Python script, data, etc. ABSOLUTE PATHS ONLY)>]```

## Advanced usage

```docker run -it --rm -p 11111:11111 -v "$HOME:/home/`whoami`/" --user $UID topologytoolkit/ttk:5.6.1-stable```

will start `pvserver` version 5.6.1 with TTK (current stable version) and listen on the default port 11111 for connections from a ParaView GUI. The directory `$(pwd)/data` will be mounted under `/home/paraview/data` in the container.

If the container is executed on a remote host, consider using the command
```
ssh -L 11111:localhost:11111 user@host docker run ...
```
which will set up the appropriate port forwarding as well. The GUI should then be able to connect to `localhost:11111`.

Notes:
- The versions of the ParaView GUI and `pvserver` have to match exactly.
- `pvserver` will currently exit after the GUI has disconnected, i.e. the container must be restarted.

## Custom Images

To re-build the image, simply clone this repository and run `docker build`. (This will use an existing docker image of ParaView and build TTK into it.)

The Dockerfile supports building a specific TTK version using the `ttk` build argument. This can be set to the designation of any branch or tag from TTK's GitHub [repository](https://github.com/topology-tool-kit/ttk), e.g. "`stable`" or "`dev`".

For example,
```
docker build -t paraview-ttk:5.8.1-dev --build-arg ttk=dev --build-arg paraview=5.8.1.
```
will build an image for TTK's `dev` branch on top of the latest ParaView release.

