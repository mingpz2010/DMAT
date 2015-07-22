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

#ifndef CROSSSECTION_H_
#define CROSSSECTION_H_

#include "Cmv.h"

// It derived from Vector_hpc, 2D ADT with { energy group + Temperature }

template<typename T> class CrossSection;
template<typename T> std::ostream& operator<<(std::ostream &s, const CrossSection<T> &M);

template<typename T>
class CrossSection : public Vector_hpc<T>
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t dim4_;
    integer_t dim5_;
    integer_t w1, w2, w3, w4;
public:
    CrossSection() : Vector_hpc<T>() {
        dim1_= dim2_ = dim3_ = dim4_ = dim5_ = 0;
        w1 = w2 = w3 = w4 = 0;
    }
    CrossSection(integer_t, integer_t, integer_t, integer_t, integer_t);
    inline const T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3, integer_t dim4, integer_t dim5) const {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3*w3+dim4*w4+dim5];
    }
    inline T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3, integer_t dim4, integer_t dim5) {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3*w3+dim4*w4+dim5];
    }
    inline T *ptr() { return Vector_hpc<T>::p_; }
    inline integer_t dim1() { return dim1_; }
    inline integer_t dim2() { return dim2_; }
    inline integer_t dim3() { return dim3_; }
    inline integer_t dim4() { return dim4_; }
    inline integer_t dim5() { return dim5_; }
    inline integer_t dim() { return Vector_hpc<T>::dim_; }
    inline void mul(T num) {
        for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
            Vector_hpc<T>::p_[i] *= num;
        }
    }

    CrossSection<T> & newsize(integer_t, integer_t, integer_t, integer_t, integer_t);

    CrossSection<T> & operator=(const CrossSection<T>&);
};

#include "CrossSection.tpp"

#endif /* CROSSSECTION_H_ */

