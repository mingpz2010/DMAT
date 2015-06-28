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
Stencil2D<T>::Stencil2D()
{
    my_func = NULL;
}

template <typename T>
Stencil2D<T>::Stencil2D(integer_t n, integer_t k)
{
    m.newsize(n, k);
    my_func = NULL;
}

template <typename T>
Stencil2D<T>::Stencil2D(stencil2d_func f, integer_t n, integer_t k)
{
    m.newsize(n, k);
    my_func = f;
}

template <typename T>
void Stencil2D<T>::stencil_reg(stencil2d_func f)
{
    my_func = f;
}

template <typename T>
void Stencil2D<T>::boundary(stencil2d_func f)
{
    if (m.dim() != 0) {
        f(m);
    }
}

template <typename T>
void Stencil2D<T>::run()
{
    if (my_func == NULL) {
        std::cerr << "Error: Stencil2D run has no stencil regedit! " << std::endl;
        return;
    }

    if (m.dim() != 0) {
        my_func(m);
    }
}



