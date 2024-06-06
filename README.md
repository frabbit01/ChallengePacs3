# ChallengePacs3
Third Challenge for the Apsc course of the 2023-2024 accademic year

Warning: 
The user can modify the json file data.json in order to change the source term for the pde or the value of n. However, if the function needs to use the value of M_PI, please write p in its place in order to get read correctly.

To make the dynamic library muParser work properly, you might need to update the path in the varialble LD_LIBRARY_PATH.

Here I report the main contents of the project files.

data.json: File that is read to initialize the forcing term, the exact solution, the number of points on each dimension of the matrix and
            Dirichlet boundary conditions (1 sx, 2 up, 3 dx, 4 down).
main.cpp: Implementation of the parallel jacobi method 
Parallel.hpp: Declaration of auxiliary functions called in the main file:
            - a local stop criterion
            - one parallel Jacobi iteration
            - L2 error
            - a function that imposes boundary condition on a matrix
Parallel.cpp: Definition of these funcions
-muparser_fun.hpp: File that contains a method that translates a string into a function
-writeVTK: File that contains methods to translate a matrix into a VTK file.
        (it was modified from the approach seen in the laboratories in order to accept Matrix objects)

-Doxyfile: File that contains indications to generate documentation correctly (used in the makefile, make doc).
The Doxygen folder contains the documentation.