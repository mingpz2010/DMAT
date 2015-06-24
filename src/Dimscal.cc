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

#include "Dimscal.h"

Dimensional_scal::Dimensional_scal(integer_t dim1, integer_t dim2, integer_t dim3)
                    : Vector_double(dim1*dim2*dim3)
{
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim_ = dim1*dim2*dim3;
    w1 = dim2*dim3;
    w2 = dim3;

    p_ = new double[dim_];
}

double Dimensional_scal::maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3)
{
    if (dim_ <= 0) {
        dim1 = dim2 = dim3 = -1;
        return 0.;
    }
    if (dim_ == 1) {
        dim1 = dim2 = dim3 = 0;
        return p_[0];
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
                if (ans < p_[index] ) {
                    ans = p_[index];
                    dim1 = i; dim2 = j; dim3 = k;
                }
            }
        }
    }

    return ans;
}

void Dimscal::copyFortran(int ref, T *, INTEGER dim1, INTEGER dim2, INTEGER dim3)
{

}




