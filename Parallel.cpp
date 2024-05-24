#include "Parallel.hpp"

namespace apsc{
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
    int 
    stop_criterion(RowMatrix &U_old, RowMatrix &U_new,const double & h,const double &tol){
        double e=0;
        auto n_r_loc=U_old.rows();
        auto n_c_loc=U_old.cols();
        for(std::size_t i=0;i<n_r_loc;++i){
            for(std::size_t j=0;j<n_c_loc;++j){
                e+=(U_new(i,j)-U_old(i,j))*(U_new(i,j)-U_old(i,j));
            }
        }
        
        e=sqrt(e*h);
        if(e<tol)
            return 1;
        return 0;

    }

    RowMatrix & parallel_jacobi(RowMatrix const & U_old,RowMatrix &U_new,const double &h, std::function<double(double,double)> const &f,const int & rank) {
        auto local_rows=U_old.rows();
        auto local_cols=U_old.cols();
        #pragma omp parallel for collapse(2)
        for(std::size_t i=1;i<rank+local_rows;++i){ //No: NEED TO MODIFY I AND J in f!
            for(std::size_t j= 1;j<local_cols;++j){
                U_new(i,j)=0.25*(U_old(i-1,j)+U_old(i+1,j)+U_old(i,j-1)+U_old(i,j+1)+h*h*f(i+rank*local_rows,j));
            }
        }
        return U_new;
    }
}