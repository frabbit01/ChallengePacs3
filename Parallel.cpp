/**
 * @file Parallel.cpp
 * @author Francesca Visalli
 * @brief Definition of the functions in Parallel.hpp.
 * @version 0.1
 * @date 2024-06-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "Parallel.hpp"
using json = nlohmann::json;
namespace apsc{
    using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
    /**
     * @brief Local stop criterion
     * 
     * @param U_old Matrix computed in the previous iteration
     * @param U_new Matrix computed in the current iteration
     * @param h Mesh spacing (assuming uniform mesh along x and y axes)
     * @param tol Tolerance
     * @return int 1 if the method has locally converged, 0 otherwise
     */
    int 
    stop_criterion(RowMatrix &U_old, RowMatrix &U_new,const double & h,const double &tol){
        double e=0;
        auto n_r_loc=U_old.rows();
        auto n_c_loc=U_old.cols();
        #pragma omp parallel for collapse(2) reduction(+:e)
        for(std::size_t i=0;i<n_r_loc;++i){
            for(std::size_t j=0;j<n_c_loc;++j){
                e+=(U_new(i,j)-U_old(i,j))*(U_new(i,j)-U_old(i,j));
            }
        }
       
        e=sqrt(e*h);
        return (e < tol) ? 1 : 0;

    }
    /**
     * @brief Computes one Jacobi parallel iteration
     * 
     * @param U_old Matrix computed in the previous iteration
     * @param U_new Matrix computed in the current iteration
     * @param h Mesh spacing (assuming uniform mesh along x and y axes)
     * @param f Forcing term
     * @param rank Rank number
     * @param size Number of ranks
     * @return RowMatrix& U_new
     */
    RowMatrix & parallel_jacobi(RowMatrix const & U_old,RowMatrix &U_new,const double &h, std::function<double(double,double)> const &f,const int & rank,const int &size) {
        auto local_rows=U_old.rows();
        auto local_cols=U_old.cols();
        RowMatrix top_row(1,local_cols);
        RowMatrix bottom_row(1,local_cols);
        int tag = 0;
        #pragma omp parallel for
        for(std::size_t j=0;j<local_cols;++j){
            top_row(0,j)=U_old(0,j);
            bottom_row(0,j)=U_old(local_rows-1,j);
        }
        if (rank < size - 1) {
        // Send last row to the next process
        MPI_Send(U_old.data() + (local_rows - 1) * local_cols, local_cols, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD);
        // Receive first row from the next process
        MPI_Recv(bottom_row.data(), local_cols, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank > 0) {
            // Receive last row from the previous process
            MPI_Recv(top_row.data(), local_cols, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Send first row to the previous process
            MPI_Send(U_old.data(), local_cols, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD);
        }

        //handle first row if not the first process
        if(rank>0){
        #pragma omp parallel for
        for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(0,j)=0.25*(top_row(0,j)+U_old(1,j)+U_old(0,j-1)+U_old(0,j+1)+h*h*f(rank*(local_rows-1)*h,j*h));
            }
        }

        //handle last row if not last process: controlla ultimo indice
        if(rank<size-1){
        #pragma omp parallel for 
        for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(local_rows-1,j)=0.25*(U_old(local_rows-2,j)+bottom_row(0,j)+U_old(local_rows-1,j-1)+U_old(local_rows-1,j+1)+h*h*f(((rank+1)*(local_rows-1)-1)*h,j*h));
            }
        }
        
        //handle internal points
        #pragma omp parallel for collapse(2)
        for(std::size_t i=1;i<local_rows-1;++i){
            for(std::size_t j= 1;j<local_cols-1;++j){
                U_new(i,j)=0.25*(U_old(i-1,j)+U_old(i+1,j)+U_old(i,j-1)+U_old(i,j+1)+h*h*f((i+rank*(local_rows-1))*h,j*h));
            }
        }
        
        return U_new;
    }
    /**
     * @brief Function that computes the L2 error and fills in the matrix corresponding to the exact solution
     * 
     * @param Global_m Matrix orresponding to the approximated solution
     * @param exact_m Matrix corresponding to the exact solution
     * @param h Mesh spacing (assuming uniform mesh along x and y axes)
     * @param n Number of points (along one of the two axes)
     * @param u_ex Exact solution
     * @return double error
     */
    double L2_error(RowMatrix & Global_m, RowMatrix & exact_m,double &h,int &n,std::function<double(double,double)> &u_ex){
        double error=0.0;
        for(int i=1;i<n-1;++i){
            for(int j=1;j<n-1;++j){
            exact_m(i,j)=u_ex(i*h,j*h); //fill the matrix
            error+=(Global_m(i,j)-exact_m(i,j))*(Global_m(i,j)-exact_m(i,j)); //update error
            }
        }
        return error=sqrt(h*error);
    }

    /**
     * @brief Fills in a matrix with four different boundary conditions taken from a json file. (1 sx, 2 up, 3 dx, 4 down)
     * 
     * @param M Matrix to be filled
     * @param h Mesh spacing (assuming a uniform mesh)
     * @return RowMatrix& M 
     */
    RowMatrix & Dirichlet_boundary_conditions(RowMatrix & M,double &h){
        std::function<double(double,double)> b_1,b_2,b_3,b_4;
        std::string string_1,string_2,string_3,string_4;
        //Read strings from json file
        std::ifstream file("data.json");
        json data = json::parse(file);
        string_1 = data.value("b_1","");
        string_2 = data.value("b_2","");
        string_3 = data.value("b_3","");
        string_4 = data.value("b_4","");
        //Initialize function
        b_1=createMuParserFunction(string_1);
        b_2=createMuParserFunction(string_2);
        b_3=createMuParserFunction(string_3);
        b_4=createMuParserFunction(string_4);

        //fill in the matrix with the boundary conditions
        for(std::size_t i=0;i<M.rows();++i){
            M(i,0)=b_1(i*h,0);
            M(i,M.cols()-1)=b_3(i*h,(M.cols()-1)*h);
        }

        for(std::size_t j=0;j<M.cols();++j){
            M(0,j)=b_2(0,j*h);
            M(M.rows()-1,j)=b_4((M.rows()-1)*h,j*h);
        }

        return M;

    }
        
}
