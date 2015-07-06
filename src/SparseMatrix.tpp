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
void SparseMatrix<T>::init(sparsematrix_manner_t type, const Matrix_hpc<T>& M)
{
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    nonzeroes = 0;

    if (i_max <=0 || j_max <= 0) {
        std::cerr << "Error: bad value in SparseMatrix constructor(<0) " << std::endl;
        exit(-1);
    }

    for (integer_t i=0; i<i_max; i++) {
        for (integer_t j=0; j<j_max; j++) {
            if (M(i,j) > 0 || M(i,j) < 0) {
                nonzeroes += 1;
            }
        }
    }

    val = NULL;
    left_val = right_val = NULL;
    col_ind = row_ptr = row_ind = col_ptr = NULL;
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
        row_ptr[i_max] = nonzeroes + 1;
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
        col_ptr[j_max] = nonzeroes + 1;
    } else if (type == TDS_MANNER) {
        this->type = 2;
        if (i_max != j_max) {
            std::cerr << "Error: bad value in SparseMatrix TDS constructor " << std::endl;
            std::cerr << " Matrix must be square matrix! "<< std::endl;
            exit(-1);
        }
        val = new T[i_max];
        left_val = new T[i_max];
        right_val = new T[i_max];

        if (i_max == 1) {
            val[0] = M(0,0);
        } else if (i_max == 2) {
            val[0] = M(0,0); val[1] = M(1,1);
            left_val[0] = 0; left_val[1] = 0;
            right_val[0] = M(0,1); right_val[1] = 0;
        } else {
            for (integer_t i=1; i<i_max-1; i++) {
                val[i] = M(i,i);
                left_val[i] = M(i, i-1);
                right_val[i] = M(i, i+1);
            }
            val[0] = M(0,0);  val[i_max-1] = M(i_max-1, i_max-1);
            left_val[0] = 0;  left_val[i_max-1] = M(i_max-1, i_max-2);
            right_val[0] = M(0, 1);  right_val[i_max-1] = 0;
        }
    } else {
        this->type = -1;
    }
}

template <typename T>
SparseMatrix<T>::SparseMatrix(Matrix_hpc<T>& M)
{
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();
    // default is CSR manner
    type = 0;
    nonzeroes = 0;
    dim1_ = i_max; dim2_ = j_max;

    if (i_max <=0 || j_max <= 0) {
        std::cerr << "Error: bad value in SparseMatrix constructor(<0) " << std::endl;
        exit(-1);
    }

    for (integer_t i=0; i<i_max; i++) {
        for (integer_t j=0; j<j_max; j++) {
            if (M(i,j) > 0 || M(i,j) < 0) {
                nonzeroes += 1;
            }
        }
    }

    val = new T[nonzeroes];
    left_val = right_val = NULL;
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
    row_ptr[i_max] = nonzeroes + 1;

    row_ind = col_ptr = NULL;
    left_val = right_val = NULL;
}

template <typename T>
SparseMatrix<T>::SparseMatrix(sparsematrix_manner_t type, Matrix_hpc<T>& M)
{
    this->init(type, M);  // M ---> const M, read-only
}

template <typename T>
void SparseMatrix<T>::chase_method(integer_t N, T *x, T *b)
{
    if (type != 2) {
        std::cerr << "Error: SparseMatrix is not TDS manner for chase method!" << std::endl;
        exit(-1);
    }

    T ans = 1.;
    for (int i=0; i<N; i++) {
        ans *= val[i];
    }
    if (ans == 0.) {
        std::cerr << "Error: SparseMatrix chase method coeff error!" << std::endl;
        exit(-1);
    }

    if (N == 1) {
        x[0] = b[0]/val[0];
        return;
    }
    if (N == 2) {
        x[1] = b[1]/val[1];
        x[0] = (b[0]-right_val[0]*x[1])/val[0];
        return;
    }

    T *beta;
    beta = new T[N];
    T *d;
    d = new T[N];
    if (beta == NULL || d == NULL) {
        std::cerr << "Error: SparseMatrix chase method memory alloc failure!" << std::endl;
        exit(-1);
    }

    beta[0] = val[0];
    x[0] = b[0];
    for (int i=1; i<N; i++) {
        d[i] = left_val[i]/beta[i-1];
        beta[i] = val[i] - d[i]*right_val[i-1];
        x[i] = b[i] - d[i]*x[i-1];
    }
    x[N-1] = x[N-1] / beta[N-1];
    for (int i=N-2; i>=0; i--) {
        x[i] = (x[i]-right_val[i]*x[i+1])/beta[i];
    }

    delete[] beta;
    delete[] d;
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const Matrix_hpc<T>& M)
{
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    type = 0;
    dim1_ = i_max;
    dim2_ = j_max;

    if (val != NULL) {
        delete[] val;
    }
    if (col_ind != NULL) {
        delete[] col_ind;
    }
    if (row_ptr != NULL) {
        delete[] row_ptr;
    }
    if (row_ind != NULL) {
        delete[] row_ind;
    }
    if (col_ptr != NULL) {
        delete[] col_ptr;
    }
    if (left_val != NULL) {
        delete[] left_val;
    }
    if (right_val != NULL) {
        delete[] right_val;
    }

    this->init(CRS_MANNER, M);

    return *this;
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
    s << "    Row Index  ";
    s.width(10);
    s << "    Col Index  ";
    s.width(20);
    s << "    Value";
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

    if (M.type_of_sparse() == 2) {
        for (integer_t k=0; k< i_max; k++) {
            s.width(10); s << k+1 <<"    ";
            s.width(10); s << k+1 <<"    ";
            s.width(20); s << M.val[k]; s << std::endl;
            if (k>0) {
                s.width(10); s << k+1 <<"    ";
                s.width(10); s << k <<"    ";
                s.width(20); s << M.left_val[k]; s << std::endl;
            }
            if (k<i_max-1) {
                s.width(10); s << k+1 <<"    ";
                s.width(10); s << k+2 <<"    ";
                s.width(20); s << M.right_val[k]; s << std::endl;
            }
        }
    }

    s << std::endl;

    return s;
}


