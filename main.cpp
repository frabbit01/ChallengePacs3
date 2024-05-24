#include "Parallel.hpp"
#include <json.hpp>
//#include <muParser.h>
#include<string>
//#include "muparser_fun.hpp"
#include<fstream>
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
    std::string funString = data.value("f","");
    //std::function<double(double,double)> f=MuparserFun(funString); //Rmk: I need to pass pi correctly!
    std::function<double(double,double)> f=[] (double x,double y) {return 8*3.14*3.14*sin(2*3.14*x)*sin(2*3.14*y);};
    int n=data.value("n",11);
    RowMatrix Global_m(n,n);
    const double h=1/(n-1);
    
    int global_convergence=0;
    tic();
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int local_nrows = (n % size > rank) ? n/size +1 : n/size; //controlla
    RowMatrix U_old(local_nrows,n), U_new(local_nrows,n);
    int local_convergence;
    //check whether parallelizing this with omp might be a good idea
    for(std::size_t k=0;k<maxit;++k){
        U_new=parallel_jacobi(U_old,U_new,h,f,rank);
        local_convergence=stop_criterion(U_old,U_new,h,tol);
        MPI_Allreduce(&local_convergence, &global_convergence, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); //I check whether the iteration has converged
        if(global_convergence==1) //If every thread converged: break
          break;
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
    MPI_Gatherv(U_new.data(), local_nrows * n, MPI_DOUBLE,
                Global_m.data(), recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    std::cout<<U_new<<std::endl;
    MPI_Finalize();
    
    if(global_convergence==1&&rank==0)
      std::cout<<"The method has converged in less than "<<maxit<<" iterations"<<std::endl;
    if(rank==0){
      std::cout<<"The resulting matrix is:"<<std::endl;
      std::cout<<Global_m<<std::endl;
    }
    return 0;
}