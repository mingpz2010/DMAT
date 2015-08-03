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
Matrix_normal<T>::Matrix_normal()
{
    type = 0;
    Vector_hpc<T>::p_ = NULL;
    Vector_hpc<T>::dim_ = 0;
    dim1_ = dim2_ = 0;
    nonzeroes = 0;
}

template <typename T>
Matrix_normal<T>::Matrix_normal(integer_t m, integer_t n)
{
    if (m<=0 || n<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = 0;
        std::cerr << "Error: bad value in Matrix_normal constructor " << std::endl;
        return;
    }
    type = 0;
    dim1_ = m;
    dim2_ = n;
    Vector_hpc<T>::p_ = new T[m*n];
    Vector_hpc<T>::dim_ = m*n;
    zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_normal<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_normal<T>::Matrix_normal(matrix_manner_t type, integer_t m, integer_t n)
{
    if (m<=0 || n<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = 0;
        std::cerr << "Error: bad value in Matrix_normal constructor " << std::endl;
        return;
    }
    if (type == ROW_MAJOR) {
        this->type = 0;
    } else {
        this->type = 1;
    }
    dim1_ = m;
    dim2_ = n;
    Vector_hpc<T>::p_ = new T[m*n];
    Vector_hpc<T>::dim_ = m*n;
    zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_normal<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_normal<T>::~Matrix_normal()
{

}

