/**
 * @file Parallel.hpp
 * @author Francesca Visalli
 * @brief Declaration of auxiliary functions for the Jacobi parallel method
 * @version 0.1
 * @date 2024-06-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<omp.h>
#include<Matrix.hpp>
#include<functional>
#include <json.hpp>
#include "muparser_fun.hpp"
#include<fstream>
namespace apsc{
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>; //I only want to split by row
    int stop_criterion(RowMatrix &U_old, RowMatrix &U_new,const double &h,const double & tol);
    RowMatrix & parallel_jacobi(RowMatrix const & U_old,RowMatrix &U_new,const double &h, std::function<double(double,double)> const &f,const int & rank,const int &size );
    double L2_error(RowMatrix & Global_m, RowMatrix & exact_m,double &h,int &n,std::function<double(double,double)> &u_ex);
    RowMatrix & Dirichlet_boundary_conditions(RowMatrix & M,double &h);
}