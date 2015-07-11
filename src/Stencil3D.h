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

#include "Dimscal.h"

template <typename T>
typedef void (*stencil3d_func)(Dimscal<T>&);

template <typename T>
class Stencil3D
{
protected:
    Dimscal<T> m;
    stencil3d_func my_func;
public:
    Stencil3D();
    Stencil3D(integer_t, integer_t, integer_t);
    Stencil3D(stencil3d_func, integer_t, integer_t, integer_t);

    void stencil_reg(stencil3d_func f);
    void stencil_reg(integer_t num, integer_t position[][3]);
    void boundary(stencil3d_func);

    void run();
};

#include "Stencil3D.tpp"

#endif /* STENCIL3D_H_ */
