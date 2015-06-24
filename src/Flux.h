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

class Flux : public Vector_hpc
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t w1, w2;
public:
    Flux() : Vector_hpc() { dim1_ = dim2_ = dim3_ = 0; w1 = w2 = 0; }
    Flux(integer_t dim1, integer_t dim2, integer_t dim3);
    inline const double& Flux::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) const {
        return p_[dim1*w1+dim2*w2+dim3];
    }
    inline double& Flux::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) {
        return p_[dim1*w1+dim2*w2+dim3];
    }
    ~Face_current();
};



#endif /* FLUX_H_ */
