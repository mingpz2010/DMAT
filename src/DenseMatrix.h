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

#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include "Matrix_hpc.h"

template <typename T>
class DenseMatrix : public Matrix_hpc<T>
{
public:
    DenseMatrix();
    DenseMatrix(integer_t, integer_t);
    DenseMatrix(matrix_manner_t, integer_t, integer_t);
};

#include "DenseMatrix.tpp"

#endif /* DENSEMATRIX_H_ */
