# Thermodynamics-Exercise--AZ420
The code used for the thermodynamics exercise. The prerequisite libraries are Numpy and open-cv.

Each file performs the following tasks:
solver_1.m:
-Performs all the calculations specified in the prompt (labelled. Comment out the parts that are not needed).
-Plots the calculated values.
-Stores the constants and parameters used by the other files.

azeotrope_finder.m:
-Finds the azeotrope composition and temperature or composition and pressure depending on settings specified in solver_1.m

cubic_eos_pure.m:
-Pure substance property calculator written by Dr. Patrick Barrie, translated into MATLAB.

dataholder.m:
-Stores the data specified in the exercise prompt.

dataholder.py:
-Stores the data specified in the exercise prompt.

RMSE_plotter.py:
-Plots the scaled value of RMSE for a range of Wilson interaction parameters.

solve_antoine.m:
-Solves the Antoine equation to find the saturation pressure or temperature.

solve_P.m:
-Fits the Wilson parameters to experimental data relating liquid composition to pressure at a constant temperature.

solve_T.m:
-Finds the value of temperature at different liquid compositions using equation 1 or 2.

wilson_temp_fitter.m:
-Fits the Wilson parameters to experimental data relating liquid composition and temperature to vapour composition at a constant pressure.
