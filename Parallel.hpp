#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<omp.h>
#include<PMatrix.hpp>

namespace apsc{
    static unsigned int maxit=1000;
    static double tol=1e-6;
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>; //I only want to split by row
    int stop_criterion(RowMatrix &U_old, RowMatrix &U_new,const double &h);
    RowMatrix & parallel_jacobi(RowMatrix const & U_old,RowMatrix &U_new,const double &h, std::function<double(double,double)> const &f );
}