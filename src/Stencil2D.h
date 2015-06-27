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

#ifndef STENCIL2D_H_
#define STENCIL2D_H_

#include "Matrix_hpc.h"

template <typename T>
typedef void (*stencil2d_func)(Matrix_hpc<T>&, integer_t, integer_t);

template <typename T>
class Stencil2D : public Matrix_hpc<T>
{
protected:
    stencil2d_func my_func;
public:
    void stencil_reg(stencil2d_func f);
    void stencil_boundary();
};

#endif /* STENCIL2D_H_ */
