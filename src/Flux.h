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

#ifndef FLUX_H_
#define FLUX_H_

#include "Cmv.h"

// It derived from Vector_hpc, 3D ADT, used for representing neutron flux

template<typename T> class Flux;
template<typename T> std::ostream& operator<<(std::ostream &s, const Flux<T> &M);

template <typename T>
class Flux : public Vector_hpc<T>
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t w1, w2;
public:
    Flux() : Vector_hpc() { dim1_ = dim2_ = dim3_ = 0; w1 = w2 = 0; }
    Flux(integer_t, integer_t, integer_t);
    inline const T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) const {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }
    inline T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }

    inline integer_t dim1() const { return dim1_; }
    inline integer_t dim2() const { return dim2_; }
    inline integer_t dim3() const { return dim3_; }

    Flux<T> & newsize(integer_t, integer_t, integer_t);
    Flux<T> & operator=(const Flux<T>&);
    Flux<T> & operator=(const T&);

    friend std::ostream& operator<< <>(std::ostream &s, const Flux<T> &M);
};

#include "Flux.tpp"

#endif /* FLUX_H_ */
