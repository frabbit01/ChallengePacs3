#include "PmatrixJ.hpp"
#include<cmath>
using namespace apsc;
using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
template <template <typename> class PMatrix>
bool 
PMatrixJ<PMatrix>:: stop_criterion(PMatrix<RowMatrix> &U_old, PMatrix<RowMatrix>&U_new,const int & n){
    double e=0;
    auto n_r_loc=U_old.get_local_nRows();
    auto n_c_loc=U_old.get_local_nCols();
    for(std::size_t i=0;i<n_r_loc;++i){
        for(std::size_t j=0;j<n_c_loc;++j){
            e+=(U_new(i,j)-U_old(i,j))*(U_new(i,j)-U_old(i,j));
        }
    }
    double h=1/(n-1);
    e=sqrt(e*h);
    if(e<tol)
        return true;
}