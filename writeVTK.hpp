/**
 * @file writeVTK.hpp
 * 
 * @brief Modified function to generate a file VTK in order to accept the Matrix type
 * @version 0.1
 * @date 2024-06-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef WRITEVTK_HPP
#define WRITEVTK_HPP

#include <vector>
#include <string>
#include <fstream>
#include<Matrix.hpp>

using namespace apsc;
using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
/**
 * @brief Function that generates a STRUCTURED VTK file with a scalar field
 * 
 * @param filename The name of the file that will be generated (if not already existing) or overwritten
 * @param scalarField The matrix from which the file will import the data
 * @param nx Biggest row index
 * @param ny Biggest column index
 * @param hx Mesh spacing (on the x-axis) assuming a uniform grid
 * @param hy Mesh spacing (on the y-axis) assuming a uniform grid
 */
void generateVTKFile(const std::string & filename, 
                     const RowMatrix & scalarField, 
                     int nx, int ny, double hx, double hy) {

    // opens the file
    std::ofstream vtkFile(filename);

    // check if the file was opened
    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    // Write VTK header
    vtkFile <<  "# vtk DataFile Version 3.0\n";
    vtkFile << "Scalar Field Data\n";
    vtkFile << "ASCII\n";                                // file format
    

    // Write grid data
    vtkFile << "DATASET STRUCTURED_POINTS\n";                             // format of the dataset
    vtkFile << "DIMENSIONS " << nx+1 << " " << ny+1 << " " << 1 << "\n";  // number of points in each direction
    vtkFile << "ORIGIN 0 0 0\n";                                          // lower-left corner of the structured grid
    vtkFile << "SPACING" << " " << hx << " " << hy << " " << 1 << "\n";   // spacing between points in each direction
    vtkFile << "POINT_DATA " << (nx+1) * (ny+1) << "\n";                  // number of points
                                                                
    
    // Write scalar field data
    vtkFile << "SCALARS scalars double\n";               // description of the scalar field
    vtkFile << "LOOKUP_TABLE default\n";                 // color table

    // Write vector field data
    for (int j = 0; j < ny+1; j++) {
        for (int i = 0; i < nx+1; i++) {
            vtkFile <<  scalarField(i,j) << "\n";
        }
    }

}

#endif // WRITEVTK_HPP