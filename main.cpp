/**
 * @file main.cpp
 * @author Francesca Visalli
 * @brief main file that performs several iterations of the Jacobi method in a parallel way and writes a VTK file from the obtained results. 
 * The forcing term,exact solution, number of nodes (equal for both x and y directions) and boundary conditions are given by the user in a json file.
 * @version 0.1
 * @date 2024-06-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "Parallel.hpp"
#include "writeVTK.hpp"
#include<string>
#include "muparser_fun.hpp"
static double c_start, c_diff;
#define tic() c_start = MPI_Wtime();
#define toc(x)                                       \
  {                                                  \
    c_diff = MPI_Wtime() - c_start;                  \
    std::cout << x << c_diff << " [s]" << std::endl; \
  }
using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
using namespace apsc;
using json = nlohmann::json;
int main(int argc, char** argv){
    int provided;
    //I initialize both n and f from a json file, so that the user can modify them as wished. 
    
    std::ifstream file("data.json");
    json data = json::parse(file);
    int n=data.value("n",11);
    std::string funString = data.value("f","");
    std::string funString2 = data.value("u_ex","");
    std::function<double(double,double)> f=createMuParserFunction(funString);
    std::function<double(double,double)> u_ex=createMuParserFunction(funString2);
    
    RowMatrix Global_m(n,n),exact_m(n,n);
    double h=1/static_cast<double>(n - 1); //mesh spacing
    unsigned iter=0; //iteration counter
    double tol=1e-12,error=0.0; //I set the tolerance
    int global_convergence=0; //flag for global convergence
    //Set boundary conditions
    Dirichlet_boundary_conditions(Global_m,h);
    Dirichlet_boundary_conditions(exact_m,h);
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int local_nrows = (n % size > rank) ? n/size +1 : n/size; 
    RowMatrix U_old(local_nrows,n), U_new(local_nrows,n);
    
    int local_convergence;
    unsigned int maxit=10000; //maximum number of iterations

    tic(); //start counting
    //Jacobi
    for(std::size_t k=0;k<maxit;++k){ 
        
        U_new=parallel_jacobi(U_old,U_new,h,f,rank,size);
        
        local_convergence=stop_criterion(U_old,U_new,h,tol); //local stop criterion
        // Synchronize all MPI processes
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&local_convergence, &global_convergence, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); //I check whether the iteration has converged
        if(global_convergence==1){
          iter=k;
          //If every thread converged: break
          break;}
        U_old=U_new; //update U_old
    }
    int* recvcounts = new int[size];
    int* displs = new int[size];
    int offset = 0;
    for (int i = 0; i < size; ++i) {
        recvcounts[i] = ((n % size > i) ? n / size + 1 : n / size) * n;
        displs[i] = offset;
        offset += recvcounts[i];
    }
    if ( rank == 0){
    toc(" matrix-vector products - Time elapsed on rank " +
          std::to_string(rank) + ": ");
  }
  //Gather all results on rank 0
    MPI_Gatherv(U_new.data(), local_nrows * n, MPI_DOUBLE,
                Global_m.data(), recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD); 
    MPI_Finalize();
    
    if(global_convergence==1&&rank==0)
      std::cout<<"The method has converged in "<<iter<<" iterations"<<std::endl;
    if(global_convergence==0&&rank==0)
      std::cout<<"The method has not converged in less than "<<maxit<<"iterations"<<std::endl;
    
    //Compute the error and fill in the matrix corresponding to the exact solution
    
    if(rank==0){
      
      error=L2_error(Global_m,exact_m,h,n,u_ex);
      std::cout<<"with a L2 error of "<<error<<std::endl;
      std::cout<<"The resulting matrix is:"<<std::endl;
      std::cout<<Global_m<<std::endl;
    
      std::cout<<"exact matrix"<<exact_m<<std::endl;

      //Generate VTKFile for both exact and approximated solutions
      generateVTKFile("mesh/approximated_laplacian.vtk", Global_m, n-1,n-1, h, h); 
      generateVTKFile("mesh/exact_laplacian.vtk", exact_m, n-1,n-1, h, h); 
    }

    return 0;
}