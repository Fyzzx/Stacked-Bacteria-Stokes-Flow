# Stacked-Bacteria-Stokes-Flow
The main file to run is titled "StackedCellsFluidForce_IB".
All other MATLAB files included here are called from this main file. 
The physics surrounding what is coded is contained in Peskin and Lim's work on immersed boundary method implementation for filaments. 

This code uses a foundation of the Immersed Boundary (IB) method stokes flow passing by two stacked cylindrical cells adhered to a surface with no-slip boundary conditions. 
The purpose of the code is to be able to allow a stokes flow to establish around two cells of user defined length, radius, and angle orientations wrt to the z-axis. 
Once the flow has been established, the user can extract the amount of force each cell exerts on the surrounding fluid and the total volume entering/exiting the eulerian grid volume
This is a foundation that can be expanded upon to include turbulence, intrinisc curvature/twist of the cells, more cells in any arrangement, and a maximum cell adhesion force prior to rupture
Can be used to extract total force on a cell to estimate cell-cell adhesin rupture strength in a simple fluid flow assay

Image of an equilibrated stokes flow around two stacked cells at user-defined angles wrt to the fluid flow. The wall or glass slide (grey) is immobile while the green cell is adhered to the wall and the red cell adhered to the green. The fluid flow vectors along one plane is shown in black.

![image](https://github.com/Fyzzx/Stacked-Bacteria-Stokes-Flow/assets/103218124/e24e1a68-d812-4814-8feb-d196d04adb96)
