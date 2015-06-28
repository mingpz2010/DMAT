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

template <typename T>
SparseMatrix<T>::SparseMatrix(Matrix_hpc<T>& M)
{
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();
    // default is CSR manner
    type = 0;
    nonzeroes = 0;
    dim1_ = i_max; dim2_ = j_max;
    for (integer_t i=0; i<i_max; i++) {
        for (integer_t j=0; j<j_max; j++) {
            if (M(i,j) > 0 || M(i,j) < 0) {
                nonzeroes += 1;
            }
        }
    }

    val = new T[nonzeroes];
    col_ind = new integer_t[nonzeroes];
    row_ptr = new integer_t[i_max + 1];

    integer_t k = 0;
    for (integer_t i=0; i<i_max; i++) {
        row_ptr[i] = k + 1;
        for (integer_t j=0; j<j_max; j++) {
            if (M(i,j) > 0 || M(i,j) < 0) {
                val[k] = M(i,j);
                col_ind[k] = j + 1;
                ++k;
            }
        }
    }
    row_ptr[i_max + 1] = nonzeroes + 1;
}

template <typename T>
SparseMatrix<T>::SparseMatrix(sparsematrix_manner_t type, Matrix_hpc<T>& M)
{
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    nonzeroes = 0;
    for (integer_t i=0; i<i_max; i++) {
        for (integer_t j=0; j<j_max; j++) {
            if (M(i,j) > 0 || M(i,j) < 0) {
                nonzeroes += 1;
            }
        }
    }

    val = col_ind = row_ptr = row_ind = col_ptr = NULL;
    dim1_ = i_max; dim2_ = j_max;

    if (type == CRS_MANNER) {
        this->type = 0;
        val = new T[nonzeroes];
        col_ind = new integer_t[nonzeroes];
        row_ptr = new integer_t[i_max + 1];

        integer_t k = 0;
        for (integer_t i=0; i<i_max; i++) {
            row_ptr[i] = k + 1;
            for (integer_t j=0; j<j_max; j++) {
                if (M(i,j) > 0 || M(i,j) < 0) {
                    val[k] = M(i,j);
                    col_ind[k] = j + 1;
                    ++k;
                }
            }
        }
        row_ptr[i_max + 1] = nonzeroes + 1;
    } else if (type == CCS_MANNER) {
        this->type = 1;
        val = new T[nonzeroes];
        row_ind = new integer_t[nonzeroes];
        col_ptr = new integer_t[j_max + 1];

        integer_t k = 0;
        for (integer_t j=0; j<j_max; j++) {
            col_ptr[j] = k + 1;
            for (integer_t i=0; i<i_max; i++) {
                if (M(i,j) > 0 || M(i,j) < 0) {
                    val[k] = M(i,j);
                    row_ind[k] = i + 1;
                    ++k;
                }
            }
        }
        col_ptr[j_max + 1] = nonzeroes + 1;
    } else if (type == BCRS_MANNER) {
        this->type = 2;
    } else {
        this->type = -1;
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const SparseMatrix<T> &M)
{
    integer_t N = M.size();
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    s << "\nNumber of Rows = " << i_max << std::endl;
    s << "\nNumber of Cols = " << j_max << std::endl;
    s << std::endl;
    s.width(10);
    s << "  Row Index ";
    s.width(10);
    s << "  Col Index ";
    s.width(20);
    s << "  Value";
    s << std::endl;

    if (M.type_of_sparse() == 0) {
        int row_index = M.p_row(0);
        for (integer_t k=0; k< N; k++) {
            if ( (k+1) >= M.p_row(row_index+1)) {
                row_index++;
            }
            s.width(10);
            s << row_index <<"    ";
            s.width(10);
            s << M.col(k) <<"    ";
            s.width(20);
            s << M(k);
            s << std::endl;;
        }
    }

    if (M.type_of_sparse() == 1) {
        int col_index = M.p_col(0);
        for (integer_t k=0; k< N; k++) {
            if ( (k+1) >= M.p_col(col_index+1)) {
                col_index++;
            }
            s.width(10);
            s << M.row(k) <<"    ";
            s.width(10);
            s << col_index <<"    ";
            s.width(20);
            s << M(k);
            s << std::endl;;
        }
    }

    s << std::endl;

    return s;
}


