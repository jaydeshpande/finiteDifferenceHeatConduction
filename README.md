# finiteDifferenceHeatConduction
Flexible finite difference transient heat conduction solver for solids

The solver provides finite difference solution to 2D conduction in solids. showPlots attribute of Heat_Equation class allows to turn real time temperature monitoring on/off.


Boundary conditions include:
1. Constant Temperature
2. Convection 
3. Mixed (Convection + Radiation) 

Usage:
1. Create Heat_Equation object with grid spacing and thermal diffusivity
2. Set boundary conditions on all 4 faces 
3. run the solver for n timesteps

# SOLVER IS NOT VALIDATED #

# Under Development #
1. Boundary conditions 
2. Exception handling 
