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
CrossSection<T>::CrossSection(integer_t m, integer_t n, integer_t k, integer_t p, integer_t q)
{
    if (m<=0 || n<=0 || k<=0 || p<=0 || q<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = 0;
        std::cerr << "Error: bad value in CrossSection constructor " << std::endl;
        return;
    }
    dim1_ = m;
    dim2_ = n;
    dim3_ = k;
    dim4_ = p;
    dim5_ = q;
    w1 = n*k*p*q;
    w2 = k*p*q;
    w3 = p*q;
    w4 = q;
    Vector_hpc<T>::p_ = new T[m*n*k*p*q];
    Vector_hpc<T>::dim_ = m*n*k*p*q;
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in CrossSection<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
}

template <typename T>
CrossSection<T>& CrossSection<T>::newsize(integer_t m, integer_t n, integer_t k,
        integer_t p, integer_t q)
{
    if (m<=0 || n<=0 || k<=0 || p<=0 || q<=0) {
        if (Vector_hpc<T>::p_ != NULL) {
            delete[] Vector_hpc<T>::p_;
        }
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = 0;
        w1 = w2 = w3 = w4 = 0;
        std::cerr << "Error: bad value in CrossSection newsize " << std::endl;
        return *this;
    }

    if (Vector_hpc<T>::p_ != NULL) {
        delete[] Vector_hpc<T>::p_;
    }
    dim1_ = m;
    dim2_ = n;
    dim3_ = k;
    dim4_ = p;
    dim5_ = q;
    w1 = n*k*p*q;
    w2 = k*p*q;
    w3 = p*q;
    w4 = q;
    Vector_hpc<T>::p_ = new T[m*n*k*p*q];
    if (Vector_hpc<T>::p_ == NULL) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = 0;
        w1 = w2 = w3 = w4 = 0;
        std::cerr << "Error: bad alloc in CrossSection newsize " << std::endl;
        return *this;
    }
    Vector_hpc<T>::dim_ = m*n*k*p*q;

    return *this;
}

template <typename T>
CrossSection<T>& CrossSection<T>::operator=(const CrossSection<T>& M)
{
    for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];
    }
}



