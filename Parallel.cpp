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

    RowMatrix & parallel_jacobi(RowMatrix const & U_old,RowMatrix &U_new,const double &h, std::function<double(double,double)> const &f,const int & rank,const int &size) {
        auto local_rows=U_old.rows();
        auto local_cols=U_old.cols();
        int last_row=0;
        RowMatrix old_row(local_cols,1);
        if(rank==size-1)
            last_row=1;
        int tag = 0;
        // Send the last row to the next process, if not the last rank
        if(rank < size - 1) {
            MPI_Send(&U_new(local_rows - 1, 0), local_cols, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD);
        }

        // Receive the last row from the previous process, if not the first rank
        if(rank > 0) {
            MPI_Recv(&old_row(0, 0), local_cols, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        #pragma omp parallel for
        for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(0,j)=0.25*(old_row(0,j)+U_old(1,j)+U_old(0,j-1)+U_old(0,j+1)+h*h*f(rank*local_rows,j));
            }
        }

        #pragma omp parallel for collapse(2)
        for(std::size_t i=1;i<local_rows-last_row;++i){
            for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(i,j)=0.25*(U_old(i-1,j)+U_old(i+1,j)+U_old(i,j-1)+U_old(i,j+1)+h*h*f(i+rank*local_rows,j));
            }
        }

        
        return U_new;
    }
}