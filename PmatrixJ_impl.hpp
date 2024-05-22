#include "PmatrixJ.hpp"

using namespace apsc;
using RowMatrix=apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR>;
template <template <typename> class PMatrix>
bool 
PMatrixJ<PMatrix>:: stop_criterion(PMatrix<RowMatrix> &U_old, PMatrix<RowMatrix>&U_new){

}