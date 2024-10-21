Coordinate grids
----------------

SIROCCO supports 3 main coordinate gridding schemes, as well as one that is
tailored to handle vertically extended disks. These schemes are

* 1-d spherical
* 2-d cylindrical
* 2-d polar or r-theta coordinates

These options are controlled by the :doc:`/input/parameters/wind/Wind/Wind.coord_system` keyword. For vertically extended disks, a modified version of the a cylindrical schme is provided where the cells viewed along the x,z plane are parallelograms, so that the grid follows the disk surface.

Although SIROCCO incorporates several models, such as the SV model for disk winds, that are continuous, the velocities and other proporites are placed on a grid, as part of the setup that goes on at the beginining of program execution.

It is up to the user to choose an appropriate coordinate system and the number of grid points 
to include in any particular run of the program.

As implemented within SIROCCO, the cells are created with a logarithmic spacing, that is the cells are larger the further they are from the central source (or disk plane).  The one exception to this is that for polar coordinates, the angular separation of cells is fixed.  For imported models, on
the other hand, the user sets the exact coordinate gridding.

Obviously, the larger the coordinate grid, 100 x 100, say, compared to 30 x 30, the better the 
grid reflects a model.  Equally obviously, the larger the coordiante grid, the larger the amount of memory the program will consume, and the larger the amount of computer time the program will take to run to completion.  The increased computing is largely associated with the fact that one needs a good representation of the spectral energy distribution in each cell in order to properly calculate the ionization state in each cell.

Although the amount of memory for particular model generaly scales with the size of the grid, different 100 x 100 models, can consume 
very different amounts of memory.  This is because for the KWD and SV parameterizations, the wind does not fill all of space.  
What really matters is the numbers of cells that are in the wind, because these are the cells for which all of the information about plasma conditions and the radition field needs to be maintained.  
So a wide angle wind with a 100 x 100 grid can take much more memory than a narrow angle wind on the same grid.

It is the number of cells that are actually in the wind that determine the fidelity of the model.  

Partial cells
===============

As note above, parameterized models often have region of space that are in the wind and regions whch are not.  If one overlays, a coordinate grid 
on such a model, there will be cells that cross edges of the wind.  These partial cells present particular problems.

In SIROCCO, velocities are interpolated on the corners of wind cells, but densities are are calculated based on 
the average radiation field in a cell, and hence 
ion densities are actually cell centered. As photons pass through a cell, they encounter resonances and the actuall opacities are 
based on an interpolated value of the densities. This presents no particular problem in regions inside the wind, but it is an issue for partial cells.

Currently, by default these cells are excluded by the calculation, and the densities of these cells are set to zero.  
Because of densities are interpolated this affects the first cell that is completely in the wind.  

There are two other alternatives:

* The partial cells can be included in the calculation.  This is reasonable for the KWD model, which has a velocity law that is easy to extend out side the wind region, but is less valid for the SV model, where this is not the case.  For the KWD model, the only issue with including partial cells is that they are "smaller" than cells which are fully in the wind, and as a result are less likely to converge as adjacent cells that are fully in the wind.

* As an advanced option, the partial cells can be excluded, but instead of setting the density of the partial cell to zero, one can assign it the density of the "nearest" cell that is totally in the wind.  In this case the ionization balance of that cell is calculated using the information in the full cell, and the ion densities that are calculated in the edge of that cell during photon transfer are just those of the center of the cell, rather than an intepolation that has on ef the endpoints at zero.  

Most of the time, the treatment of partial cells does not change the predicted spectrum significantly, but this is something that is worthwhiled chacking. Users should be wary in situations where there are directions in which significant numbers of photons will pass though very few cells in the wind.  This could happen for a "narrow" wind with a very small opening angle.  Having a small number of cells in the wind is, of course, one should be concerned about in any event.  
