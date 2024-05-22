#ifndef PMATRIXJ_HPP
#define PMATRIXJ_HPP
#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<omp.h>
#include"PMatrix_modified.hpp"
#include<functional>
namespace apsc{
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>; //I only want to split by row
   template <template <typename> class PMatrix> class PMatrixJ{ //controllare se così facendo posso accedere a tutti i membri di PMatrix
        public:
            bool stop_criterion(PMatrix<RowMatrix> &U_old, PMatrix<RowMatrix>&U_new,const int & n);
            PMatrix<RowMatrix> & parallel_jacobi(PMatrix<RowMatrix> & U0, const int & n,std::function<double(double)> const &f );

        private:
            static unsigned int maxit=1000;
            static double tol=1e-6;
    };
}

#include "PmatrixJ_impl.hpp"

#endif /*PMATRIXJ_HPP*/