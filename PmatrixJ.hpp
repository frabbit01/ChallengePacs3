#ifndef PMATRIXJ_HPP
#define PMATRIXJ_HPP
#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<omp.h>
#include<PMatrix.hpp>

namespace apsc{
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>; //I only want to split by row
   template <template <typename> class PMatrix> class PMatrixJ{ //controllare se così facendo posso accedere a tutti i membri di PMatrix
        public:
            bool stop_criterion(PMatrix<RowMatrix> &U_old, PMatrix<RowMatrix>&U_new);
            PMatrix<RowMatrix> & parallel_jacobi(PMatrix<RowMatrix> & U0, const int & maxit,const double & tol);
    };
}

#include "PmatrixJ_impl.hpp"

#endif /*PMATRIXJ_HPP*/