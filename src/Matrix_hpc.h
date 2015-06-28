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
    integer_t dim1_;
    integer_t dim2_;
    integer_t nonzeroes;
public:
    Matrix_hpc();
    Matrix_hpc(integer_t, integer_t);
    Matrix_hpc(matrix_manner_t, integer_t, integer_t);
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

    inline integer_t size() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim1() const { return dim1_;}
    inline integer_t dim2() const { return dim2_;}
    inline integer_t ref() const { return Vector_hpc<T>::ref_; }
    inline int null() const {return Vector_hpc<T>::dim_== 0;}
    inline void set_nonzero(integer_t n) { nonzeroes = n; }

    Matrix_hpc<T> & newsize(integer_t, integer_t);
    Matrix_hpc<T> & operator=(const Matrix_hpc<T>&);
    Matrix_hpc<T> & operator=(const T&);

    void add(const Matrix_hpc<T> &c1);

    friend std::ostream& operator<< <>(std::ostream &s, const Matrix_hpc<T> &M);
};

#include "Matrix_hpc.tpp"

#endif /* MATRIX_HPC_H_ */
