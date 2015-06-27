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

#ifndef DIMSCAL_H_
#define DIMSCAL_H_

#include "Cmv.h"

// It derived from Vector_hpc, 3D ADT

template<typename T> class Dimscal;
template<typename T> std::ostream& operator<<(std::ostream &s, const Dimscal<T> &M);

template <typename T>
class Dimscal : public Vector_hpc<T>
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t w1, w2;
public:
    Dimscal() : Vector_hpc<T>() { dim1_= dim2_ = dim3_ = 0; w1 = w2 = 0; }
    Dimscal(integer_t, integer_t, integer_t);
    inline const T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) const {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }
    inline T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }

    // common functions
    T maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3);

    friend std::ostream& operator<< <>(std::ostream &s, const Dimscal<T> &M);
};

#include "Dimscal.tpp"

#endif /* DIMSCAL_H_ */
