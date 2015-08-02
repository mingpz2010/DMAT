/*
    Copyright (C) <2014-2020>  <PingzhouMing>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef MATRIX_HPC_H_
#define MATRIX_HPC_H_

#include "Cmv.h"

typedef enum Matrix_store_manner {
    ROW_MAJOR,
    COL_MAJOR
}matrix_manner_t;

typedef enum Matrix_property {
    M_NORMAL,
    M_BANDED,
    M_SYMMETRIC,
    M_SYMMETRIC_BANDED,
    M_TRIANGULAR
}matrix_property_t;

// " More C++ Idioms/Making New Friends "
// From the <More C++ Idioms> wiki book
// Declare friend function
template<typename T> class Matrix_hpc;
template<typename T> std::ostream& operator<<(std::ostream &s, const Matrix_hpc<T> &M);

template <typename T>
class Matrix_hpc : public Vector_hpc<T>
{
protected:
    int type;  // 0-ROW_MAJOR, default, 1-COL_MAJOR
    int property;
    integer_t dim1_;
    integer_t dim2_;
    integer_t ksub;
    integer_t ksupper;
    integer_t nonzeroes;
public:
    Matrix_hpc();
    Matrix_hpc(integer_t, integer_t);
    Matrix_hpc(matrix_manner_t, integer_t, integer_t);
    ~Matrix_hpc();
    inline T& operator()(integer_t m, integer_t n) {
        if (type == 0) {
            return Vector_hpc<T>::p_[m*dim2_ + n];
        } else {
            return Vector_hpc<T>::p_[m + n*dim1_];
        }
    }
    inline const T& operator()(integer_t m, integer_t n) const {
        if (type == 0) {
            return Vector_hpc<T>::p_[m*dim2_ + n];
        } else {
            return Vector_hpc<T>::p_[m + n*dim1_];
        }
    }

    inline T *ptr() { return Vector_hpc<T>::p_; }
    inline integer_t size() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim1() const { return dim1_;}
    inline integer_t dim2() const { return dim2_;}
    inline integer_t ref() const { return Vector_hpc<T>::ref_; }
    inline int null() const {return Vector_hpc<T>::dim_== 0;}
    inline void zero() {
        for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
            Vector_hpc<T>::p_[i] = 0;
        }
    }
    inline void set_nonzero(integer_t n) { nonzeroes = n; }

    // BLAS cutting and implementation
    void blas_op(T a, T b, T c);
    // BLAS-1
    void dswap(const Matrix_hpc<T>& M);
    void dscal(T a);
    void dcopy(const Matrix_hpc<T>& M);
    void daxpy(T a, const Matrix_hpc<T>& M);
    T dnrm2();
    T dasum();

    // BLAS-2
    void dgemv(Vector_hpc<T>& v);     // v = M*v
    void dgemv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2);   // v2 = M*v1
    void dgbmv(Vector_hpc<T>& v);
    void dgbmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2);
    void dsymv(Vector_hpc<T>& v);
    void dsymv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2);
    void dsbmv(Vector_hpc<T>& v) { dgbmv(v); }
    void dsbmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2) { dgbmv(v1, v2); }
    void dspmv(Vector_hpc<T>& v) { dsymv(v); }
    void dspmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2) { dsymv(v1, v2); }
    void dtrmv(Vector_hpc<T>& v);
    void dtrmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2);
    void dtrsv(Vector_hpc<T>& x, const Vector_hpc<T>& b);

    // BLAS-3
    void dgemm(Matrix_hpc<T>& m);       // m = M*m
    void dgemm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2);  // m2 = M * m1
    void dsymm(Matrix_hpc<T>& m);
    void dsymm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2);
    void dtrmm(Matrix_hpc<T>& m);
    void dtrmm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2);

    Matrix_hpc<T> & newsize(integer_t, integer_t);
    Matrix_hpc<T> & operator=(const Matrix_hpc<T>&);
    Matrix_hpc<T> & operator=(const T&);

    void add(const Matrix_hpc<T> &c1);

    friend std::ostream& operator<< <>(std::ostream &s, const Matrix_hpc<T> &M);
};

#include "Matrix_hpc.tpp"

#endif /* MATRIX_HPC_H_ */
