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
LegendreScal<T>::LegendreScal(integer_t m, integer_t n, integer_t k,
        integer_t p, integer_t q, integer_t r)
{
    if (m<=0 || n<=0 || k<=0 || p<=0 || q<=0 || r<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = dim6_ = 0;
        std::cerr << "Error: bad value in LegendreScal constructor " << std::endl;
        return;
    }
    dim1_ = m;
    dim2_ = n;
    dim3_ = k;
    dim4_ = p;
    dim5_ = q;
    dim6_ = r;
    w1 = n*k*p*q*r;
    w2 = k*p*q*r;
    w3 = p*q*r;
    w4 = q*r;
    w5 = r;
    Vector_hpc<T>::p_ = new T[m*n*k*p*q*r];
    Vector_hpc<T>::dim_ = m*n*k*p*q*r;
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in LegendreScal<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
}

template <typename T>
LegendreScal<T>& LegendreScal<T>::newsize(integer_t m, integer_t n, integer_t k,
        integer_t p, integer_t q, integer_t r)
{
    if (m<=0 || n<=0 || k<=0 || p<=0 || q<=0) {
        if (Vector_hpc<T>::p_ != NULL) {
            delete[] Vector_hpc<T>::p_;
        }
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = dim6_ = 0;
        w1 = w2 = w3 = w4 = w5 = 0;
        std::cerr << "Error: bad value in LegendreScal newsize " << std::endl;
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
    dim6_ = r;
    w1 = n*k*p*q*r;
    w2 = k*p*q*r;
    w3 = p*q*r;
    w4 = q*r;
    w5 = r;
    Vector_hpc<T>::p_ = new T[m*n*k*p*q*r];
    if (Vector_hpc<T>::p_ == NULL) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = dim3_ = dim4_ = dim5_ = dim6_ = 0;
        w1 = w2 = w3 = w4 = w5 = 0;
        std::cerr << "Error: bad alloc in LegendreScal newsize " << std::endl;
        return *this;
    }
    Vector_hpc<T>::dim_ = m*n*k*p*q*r;

    return *this;
}

template <typename T>
LegendreScal<T>& LegendreScal<T>::operator=(const LegendreScal<T>& M)
{
    for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];
    }
}



