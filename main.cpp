#include "Parallel.hpp"
#include <json.hpp>
//#include <muParser.h>
#include<string>
#include "muparser_fun.hpp"
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
    //I initialize both n and f from a json file, so that the user can modify them as wished. MUPARSER NOT WORKING!!!
    std::ifstream file("data.json");
    json data = json::parse(file);
    std::string funString = data.value("f","");
    std::cout<<funString<<std::endl;
    //std::function<double(double,double)> f=MuparserFun(funString); //Rmk: I need to pass pi correctly!
    std::function<double(double,double)> f=[] (double x,double y) {return 8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);};
    int n=data.value("n",11);
    RowMatrix Global_m(n,n);
    double h=1/static_cast<double>(n - 1);
    unsigned iter=0;
    double tol=1e-7;
    int global_convergence=0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int local_nrows = (n % size > rank) ? n/size +1 : n/size; //controlla
    RowMatrix U_old(local_nrows,n), U_new(local_nrows,n);
    // Initialize U_old with a non-trivial guess
    // Create a random device and seed the random number generator
    /*std::random_device rd;
    std::mt19937 gen(rd());

    // Define the range for the random number
    std::uniform_real_distribution<> distrib(-10, 10);

    // Generate and print a random number
    
    for (std::size_t i = 0; i < local_nrows; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          int random_number = distrib(gen);
            U_old(i, j) = random_number;
        }
    }*/

    int local_convergence;
    //check whether parallelizing this with omp might be a good idea
    unsigned int maxit=10000;
    tic(); //start counting

    for(std::size_t k=0;k<maxit;++k){ //immediately gets out??? maybe 0 is a solution for this problem??
        
        U_new=parallel_jacobi(U_old,U_new,h,f,rank,size);
        
        local_convergence=stop_criterion(U_old,U_new,h,tol);
        // Synchronize all MPI processes
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&local_convergence, &global_convergence, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); //I check whether the iteration has converged
        if(global_convergence==1){
          iter=k;
           //If every thread converged: break
          break;}
        U_old=U_new; //update U_old
        //MPI_Barrier(MPI_COMM_WORLD);
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
    //std::cout<<U_new<<std::endl;
    MPI_Finalize();
    
    if(global_convergence==1&&rank==0)
      std::cout<<"The method has converged in "<<iter<<" iterations"<<std::endl;
    if(global_convergence==0&&rank==0)
      std::cout<<"The method has not converged in less than "<<maxit<<"iterations"<<std::endl;
    if(rank==0){
      std::cout<<"The resulting matrix is:"<<std::endl;
      std::cout<<Global_m<<std::endl;
    }

    
    return 0;
}