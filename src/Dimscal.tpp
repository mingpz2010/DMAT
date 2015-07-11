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
Dimscal<T>::Dimscal(integer_t m, integer_t n, integer_t k)
{
    if (m<=0 || n<=0 || k<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = 0;
        std::cerr << "Error: bad value in Dimscal constructor " << std::endl;
        return;
    }
    dim1_ = m;
    dim2_ = n;
    dim3_ = k;
    w1 = n*k;
    w2 = k;
    Vector_hpc<T>::p_ = new T[m*n*k];
    Vector_hpc<T>::dim_ = m*n*k;
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Dimscal<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
}

template <typename T>
void Dimscal<T>::blas_op(T a, T b, T c)
{
    if (c == 0)  return;

    for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
        Vector_hpc<T>::p_[i] *= a;
        Vector_hpc<T>::p_[i] += b;
        Vector_hpc<T>::p_[i] /= c;
    }
}

template <typename T>
T Dimscal<T>::maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3)
{
    if (Vector_hpc<T>::dim_ <= 0) {
        dim1 = dim2 = dim3 = -1;
        return 0.;
    }
    if (dim_ == 1) {
        dim1 = dim2 = dim3 = 0;
        return Vector_hpc<T>::p_[0];
    }

    integer_t i,j,k;
    integer_t index;
    double ans;

    ans = 0.;
    for (i=0; i<dim1_; i++) {
        index = i*w1;
        for (j=0; j<dim2_; j++) {
            index += j*w2;
            for (k=0; k<dim3_; k++) {
                index += k;
                if (ans < Vector_hpc<T>::p_[index] ) {
                    ans = Vector_hpc<T>::p_[index];
                    dim1 = i; dim2 = j; dim3 = k;
                }
            }
        }
    }

    return ans;
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const Dimscal<T> &M)
{
    integer_t N = M.size();
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();
    integer_t k_max = M.dim3();

    for (integer_t i=0; i< i_max; i++) {
        for (integer_t j=0; j< j_max; j++) {
            for (integer_t k=0; k< k_max; k++) {
                s << M(i, j, k) << std::endl;
            }
        }
    }

    s << std::endl;

    return s;
}



