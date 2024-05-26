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
        RowMatrix top_row(1,local_cols);
        RowMatrix bottom_row(1,local_cols);
        int last_row = (rank == size - 1) ? 1 : 0;
        int tag = 0;
        #pragma omp parallel for
        for(std::size_t j=0;j<local_cols;++j){
            top_row(0,j)=U_old(0,j);
            bottom_row(0,j)=U_old(local_rows-1,j);
        }
        //sistemare i send e receive
        if (rank > 0) {
        // Send the first row to the next process, if not the first rank
        MPI_Send(top_row.data(), local_cols, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD);
        // Receive the last row from the previous process, if not the first rank
        MPI_Recv(top_row.data(), local_cols, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        
        if (rank < size - 1) {
            // Send the last row to the next process, if not the last rank
            MPI_Send(bottom_row.data(), local_cols, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD);
            // Receive the first row from the next process, if not the last rank
            MPI_Recv(bottom_row.data(), local_cols, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(rank>0){
        #pragma omp parallel for
        for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(0,j)=0.25*(top_row(0,j)+U_old(1,j)+U_old(0,j-1)+U_old(0,j+1)+h*h*f(rank*local_rows,j));
            }
        }

        if(rank<size-1){
        #pragma omp parallel for
        for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(local_rows-1,j)=0.25*(U_old(local_rows-1,j)+bottom_row(0,j)+U_old(local_rows-1,j-1)+U_old(local_rows-1,j+1)+h*h*f(rank*local_rows,j));
            }
        }
        
        //handle internal points
        #pragma omp parallel for collapse(2)
        for(std::size_t i=1;i<local_rows-1;++i){
            for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(i,j)=0.25*(U_old(i-1,j)+U_old(i+1,j)+U_old(i,j-1)+U_old(i,j+1)+h*h*f(i+rank*(local_rows-1),j));
            }
        }
        
        return U_new;
    }
}
