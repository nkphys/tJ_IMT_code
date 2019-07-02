#include "Matrix.h"
#include "tensor_type.h"


#ifndef define_functions
#define define_functions

Matrix<complex<double>> product(Matrix<complex<double>> A, Matrix<complex<double>> B){
    assert(A.n_col()==B.n_row());

    Matrix<complex<double>> C;
    C.resize(A.n_row(),B.n_col());

    for(int i=0; i<A.n_row(); i++) {
        for(int j=0; j<B.n_col(); j++) {

            C(i,j) = zero_complex;
            for(int k=0; k<B.n_row(); k++) {
            C(i,j) += A(i,k)*B(k,j);
            }

        }
    }

    return C;
} // ----------


complex<double> DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right){
    complex<double> temp_;
    temp_=zero_complex;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += conj(left[i])*right[i];
    }
    return temp_;

}

double DOT_P(Mat_1_doub left, Mat_1_doub right){
    double temp_;
    temp_=0.0;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += left[i]*right[i];
    }
    return temp_;

}


Matrix<double> product(Matrix<double> A, Matrix<double> B){
    assert(A.n_col()==B.n_row());

    Matrix<double> C;
    C.resize(A.n_row(),B.n_col());

    for(int i=0; i<A.n_row(); i++) {
        for(int j=0; j<B.n_col(); j++) {

            C(i,j) = 0.0;
            for(int k=0; k<B.n_row(); k++) {
            C(i,j) += A(i,k)*B(k,j);
            }

        }
    }

    return C;
} // ----------


#endif
