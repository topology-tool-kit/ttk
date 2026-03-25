# TTK DebrisTracer — README

## Overview

**TTK DebrisTracer** post-processes trajectories extracted from time-varying scalar fields to:
- **linearize** each trajectory segment,
- **fuse/chain** compatible trajectories into longer **merged chains** using direction, time-gap, and distance constraints,
- **extract surface regions** on a selected frame (BFS / Random-Walker) and attach them back to trajectories,
- **compute per-chain statistics** (duration, velocity, surface metrics).

This filter is intended to be used **after** TTK’s tracking filter *TrackingFromFields*, together with the original dataset carrying point scalars for each frame. 

> **Dependency:** this filter **requires Eigen** for linear regression and Random-Walker segmentation. Make sure TTK is built with `TTK_ENABLE_EIGEN=ON`.

---

## Inputs

The filter takes **two inputs** 

1) **InputGrid (port 0)** — `vtkUnstructuredGrid`  
   Output of **TrackingFromField** (or equivalent): each **cell** represents a **trajectory segment** (typically a line between two successive points). Expected arrays include:
   - Per-cell **ConnectedComponentId** (trajectory id),
   - Per-point **TimeStep** (frame index),
   - Per-point **VertexGlobalId** (vertex id in the original mesh),
   - Per-point **InstantPersistence**.

2) **InputSet (port 1)** — `vtkDataSet`  
   The **original time-slice dataset** carrying **point scalar fields** (one array per frame) used for surface extraction and gradient magnitudes.

---

## Parameters

### Fusion / chaining
Controls how independent trajectories are **linked** into longer **chains**.

- **Chaining: Min cosine similarity** (`cosCol`)  
  Minimal cosine between **mean unit directions** of two trajectories to consider linking.

- **Chaining: Max 3D distance²** (`maxRadus`)  
  Maximum allowed **squared distance** in **(x,y,t)** at the candidate linking frame between the **predicted end of i** and the **start of j**.

- **Chaining: frame gap** (`maxFrameDist`, with implicit lower bound)  
  Frame‐gap window between `end(i)` and `start(j)` for valid links : [-value,value].

### Surface extraction
Extracts, on a selected frame, **surface regions** around each trajectory point and aggregates statistics.

- **Surface: Frame index** (`frameSurface`)  
  The frame at which surfaces are displayed.

- **Segment cleanup threshold** (`errSurf`)
  Minimum number of Otsu bins for dark-segment cleanup (0 = disabled).- **Surface: Method** (`surfaceMethod`, advanced) 
 
- **Persistence Threshold (%)** (`persisThresh`)
  Persistence simplification threshold as a percentage of maximum persistence.

### Units
Unit conversion and duration reporting. 

- **Units: Spatial scale (px→mm)** (`spatialScale`)  
  Millimeters per pixel; used to convert velocities.

- **Units: Inter-frame time (μs)** (`interFrame`)  
  Time step between frames, in microseconds.

- **Units: Durations in seconds** (`convertDur`, advanced, boolean)  
  If enabled, durations are reported in **seconds** instead of frames.

### Directional filtering
Optional additional filters on direction/velocity. 

- **Filter: Min |Vx| (mm/s)** (`minVx`)  
  Minimal absolute X-velocity (after unit conversion) for a trajectory to be accepted.

- **Direction: Y cosine max** (`filtreY`)  
  Enforces a symmetric bound on the normalized Y component of direction (|cos(Y)| ≤ value).  

### Initial position filtering
Filters trajectories **based on their position at the first frame** (initial position). Use `-1` to disable a bound.

- **Initial position: Min X / Max X** (`minX`, `maxX`)  
  Lower/upper bound on **X** at the trajectory’s initial frame (Paraview coordinates).
- **Initial position: Min Y / Max Y** (`minY`, `maxY`)  
  Lower/upper bound on **Y** at the trajectory’s initial frame (Paraview coordinates).

### Advanced drawing
- **Extend merged chains to t=0** (`extendTraj`, boolean)  
  If ON, merged chains are drawn starting at **frame 0**; otherwise, from their own start frame. 

---

## Outputs

The filter produces **four outputs** 

1) **Statistics (port 0)** — `vtkTable`  
   One row per **merged chain** with the following columns:
   - **StartFrame** / **EndFrame** — chain time bounds.  
   - **Duration** — chain length (in frames or seconds if `convertDur` is ON).  
   - **VX**, **VY** — average planar velocity components.  
   - **SurfaceMean** — mean surface size aggregated over member trajectories.  
   - *(optionally present if enabled in your build)* **SurfaceMin**, **SurfaceMax** — per-chain min/max surface sizes.

2) **Trajectories (port 1)** — `vtkUnstructuredGrid`  
   - **Traj Id** — id of the final merged chain 
   - **Duration** — endFrame - startFrame

3) **Surfaces (port 2)** — `vtkDataSet` (same type as the second input)  
   Copy of the original dataset **augmented** with a point-data array:
   - **FinalTrajId** — for each vertex:  
     - `≥ 0`: id of the **final merged chain** that claimed the vertex on `frameSurface`,  
     - `-1`: vertex **not part of any** extracted surface,  
     - `-2`: vertex is a **duplicate** (claimed by multiple trajectories).

4) **Chains (port 3)** — `vtkUnstructuredGrid`  
   Contains both: the per-trajectory **linearized segments** (initial fits), the additional **“fusion links”** (segments connecting i → j).
   - **FinalChainId** : final merged chain id carried by the contributing trajectory
   - **InputTrajId**  : original trajectory id for initial segments, -1 for links
   - **SegmentKind**  : 0 = initial linear segment, 1 = fusion link




