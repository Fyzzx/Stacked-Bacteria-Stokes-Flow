# Stacked-Bacteria-Stokes-Flow
Uses a foundation of an Immersed Boundary method flow passing by two stacked cylindrical cells adhered to a surface with no-slip boundary conditions. 
The purpose of the code is to be able to allow a stokes flow to establish around two cells of user defined length, radius, and angle orientations wrt to the z-axis. 
Once the flow has been established, the user can extract the amount of force each cell exerts on the surrounding fluid and the total volume entering/exiting the eulerian grid volume
This is a foundation that can be expanded upon to include turbulence, intrinisc curvature/twist of the cells, more cells in any arrangement, and a maximum cell adhesion force prior to rupture
Can be used to extract total force on a cell to estimate cell-cell adhesin rupture strength in a simple fluid flow assay
