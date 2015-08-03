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
#ifndef MATRIX_NORMAL_H_
#define MATRIX_NORMAL_H_

#include "Common.h"
#include "Cmv.h"

template <typename T>
class Matrix_normal : public Vector_hpc<T>
{
protected:
    int type;  // 0-ROW_MAJOR, default, 1-COL_MAJOR
    integer_t dim1_;
    integer_t dim2_;
    integer_t nonzeroes;
public:
    Matrix_normal();
    Matrix_normal(integer_t, integer_t);
    Matrix_normal(matrix_manner_t, integer_t, integer_t);
    ~Matrix_normal();
    inline T *ptr() { return Vector_hpc<T>::p_; }
    inline integer_t size() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim() const { return Vector_hpc<T>::dim_;}
    inline integer_t dim1() const { return dim1_;}
    inline integer_t dim2() const { return dim2_;}
    inline int null() const {return Vector_hpc<T>::dim_== 0;}
    inline void zero() {
        for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
            Vector_hpc<T>::p_[i] = 0;
        }
    }
};


#include "Matrix_normal.tpp"

#endif /* MATRIX_NORMAL_H_ */
