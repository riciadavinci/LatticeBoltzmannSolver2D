# SIWIR2 Assignment 4

# Team Members:
1. Nikhil VC
2. Rishikesh AN
3. Nagesh AC

# Implementing Lattice Boltzmann Solver for Obstacle in a Channel

# Compilation (makes use of g++ & c++17 standard):
make

# Usage:
    ./lbm <PARAMETERS>.dat

Here:\n"
    => <PARAMETERS>.dat : parameters file


# Usage examples:
1. ./lbm files/params_Re40.dat
2. ./lbm files/params_Re100_wing.dat
3. ./lbm files/params_Re500.dat

# Implemented the bonus task as well: Re100_wing.dat

# Note:
If you want to include a similar PARAMETERS.dat file like Re100_wing.dat, where the geometry is 
passed via the file, make sure that the geometry file <GEOMETRY>.pgm is located in the "files" 
directory within the project and also the geometry parameter should just be the filename:
# PARAMETERS.dat file
...
geometry <some_geometry>.pgm
...

# Utility:
use "make vtkclean" to delete/clean all .vtk files generated within the directory
