# Paraview + TTK Docker Image

This docker file can be used to build docker images containing installations of the [Topology Tool Kit (TTK)](http://topology-tool-kit.github.io) and/or [ParaView](http://www.paraview.org) and corresponding dependencies for usage and development. For a description of the build process, see [Custom Images](##custom-images), for usageof the images, see [Simple Usage](##simple-usage) or [Advanced Usage](##advanced-usage).

The full docker image specifically contains:

- ParaView server with offscreen rendering using either [OSMesa](http://www.mesa3d.org/osmesa.html) or [OSPRay](http://www.ospray.org).
- TTK for ParaView plugins are installed.

It is supposed to be used in conjunction with a local ParaView GUI.

## Simple usage

To run Kitware's binary distribution of ParaView with TTK's docker, simply run:

```
./runParaViewTTKDocker.sh <Path to ParaView binary (version 5.10.1)> [<Standard ParaView arguments (state files, data, etc.)>]
```

To run a python script which uses TTK, simply run:

```
./runTTKPythonDocker.sh [<Standard pvpython arguments: Python script, data, etc. ABSOLUTE PATHS ONLY)>]
```

## Custom Images

To re-build the image, simply clone this repository and run `docker build .` (which will build a docker image with current versions of ParaView and TTK).

The Dockerfile can build different images using the following targets through `docker build --target <targetname> .`:
- `ttk` will build the default target containing ParaView and TTK. When run, the image starts a ParaView server on port 11111.
- `ttk-dev` will build a docker image containing an installation of ParaView built from source and development tools such as cmake. It can be used for development of your own version of TTK.

The Dockerfile supports building a specific TTK version using the `ttk` build argument. This can be set to the designation of any branch or tag from TTK's GitHub [repository](https://github.com/topology-tool-kit/ttk), e.g. "`master`" or "`v0.9.7`". The default value is the current ttk dev branch.

The Dockerfile supports building a specific ParaView version using the `paraview` build argument. This can be set to any release tag from ParaViews's [download page](https://www.paraview.org/download/), e.g. "`5.10.0`" or "`5.9.1`". Note however, that versions older than `5.9` are no longer maintained. The default value is the latest ParaView release.

For example, `docker build -t ttk --build-arg paraview=5.9.1 --build-argttk=1.0.0 .` builds a docker image which by default will run a ParaView 5.9.1 server containing a TTK 1.0.0 installation.

## Advanced usage

After building the image with `docker build -t <image-name> --build-arg paraview=<pv-version> --build-argttk=<ttk-version> .` you can run it with
```
docker run -it --rm -p 11111:11111 -v "$HOME:/home/`whoami`/" --user $UID <image-name>
```
or use the GitHub image with 
```
docker run -it --rm -p 11111:11111 -v "$(pwd)/data:/home/paraview/data" --user $UID ghcr.io/scivislab/ttk:latest
```

which will start `pvserver` (version `<pv-version>`) with TTK (version <ttk-version>) and listen on the default port 11111 for connections from a ParaView GUI. The directory `$(pwd)/data` will be mounted under `/home/paraview/data` in the container.

If the container is executed on a remote host, consider using the command
```
ssh -L 11111:localhost:11111 user@host docker run ...
```
which will set up the appropriate port forwarding as well. The GUI should then be able to connect to `localhost:11111`.

Notes:
- The versions of the ParaView GUI and `pvserver` have to match exactly.
- `pvserver` will currently exit after the GUI has disconnected, i.e. the container must be restarted.
