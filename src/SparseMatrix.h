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

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include "Matrix_hpc.h"

typedef enum Sparse_matrix_store_manner {
    CRS_MANNER,         // Compressed Row Storage
    CCS_MANNER,         // Compressed Column Storage
    TDS_MANNER,         // Three Diagonal Storage, Sparse Matrix must be square matrix
}sparsematrix_manner_t;

template<typename T> class SparseMatrix;
template<typename T> std::ostream& operator<<(std::ostream &s, const SparseMatrix<T> &M);

template <typename T>
class SparseMatrix
{
private:
    void init(sparsematrix_manner_t type, const Matrix_hpc<T>&);
protected:
    int type;  // 0-CRS, default, 1-CCS, 2-CDS
    integer_t dim1_;
    integer_t dim2_;
    integer_t nonzeroes;
    T *val;
    // CRS manner
    integer_t *col_ind;
    integer_t *row_ptr;
    // CCS manner
    integer_t *row_ind;
    integer_t *col_ptr;
    // TDS manner
    T * left_val;
    T * right_val;
public:
    SparseMatrix() {
        type = nonzeroes = 0;
        dim1_ = dim2_ = 0;
        val = left_val = right_val = NULL;
        col_ind = row_ptr = row_ind = col_ptr = NULL;
    }
    SparseMatrix(Matrix_hpc<T>&);
    SparseMatrix(sparsematrix_manner_t type, Matrix_hpc<T>&);
    inline T& operator()(integer_t k) {
        return val[k];
    }
    inline const T& operator()(integer_t k) const {
        return val[k];
    }

    inline integer_t size() const { return nonzeroes;}
    inline integer_t dim1() const { return dim1_;}
    inline integer_t dim2() const { return dim2_;}
    // CRS manner
    inline integer_t col(integer_t k) const { return col_ind[k]; }
    inline integer_t p_row(integer_t k) const { return row_ptr[k]; }
    // CCS manner
    inline integer_t row(integer_t k) const { return row_ind[k]; }
    inline integer_t p_col(integer_t k) const { return col_ptr[k]; }
    inline int type_of_sparse() const { return type; }

    // In order to unified, CRS manner is stored
    SparseMatrix<T> & operator=(const Matrix_hpc<T>&);

    friend std::ostream& operator<< <>(std::ostream &s, const SparseMatrix<T> &M);
};

#include "SparseMatrix.tpp"

#endif /* SPARSEMATRIX_H_ */
