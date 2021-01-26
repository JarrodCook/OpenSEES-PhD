## Generic linear dynamic setup for a 2D structure using Finite Element Analysis written in MATLAB

### Files:

#### generic_linear_dynamic.m
-Main file to enter settings and run the analysis, calls upon the following function files

#### Assembly.m
-Create assembly matrices from dof matches

#### anglefun.m
-Returns element angles and lengths from coordinates

#### frame.m
-Return stiffness properties for frame element

#### dynamic1.m
-Newmark-beta incremental integration algorithm

#### fordis2.m
-Returns forces and displacements in global and local coordinates for each element

#### deflection2.m
-Produce coordinates of deformed node positions from original coordinates and FEA displacements

#### eulersf3.m
-Returns element displacement between nodes using Euler beam shape functions
