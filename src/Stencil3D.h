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

#ifndef STENCIL3D_H_
#define STENCIL3D_H_

#include "Cmv.h"

template <typename T>
typedef void (*stencil3d_func)(T ***, integer_t, integer_t, integer_t);

template <typename T>
class Stencil3D : public Vector_hpc<T>
{
protected:
    stencil3d_func my_func;
public:
    void stencil_reg(stencil3d_func f);
    void stencil_boundary();
};


#endif /* STENCIL3D_H_ */
