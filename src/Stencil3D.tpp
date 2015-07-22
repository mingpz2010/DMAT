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
Stencil3D<T>::Stencil3D()
{
    my_func = NULL;
}

template <typename T>
Stencil3D<T>::Stencil3D(integer_t n, integer_t k, integer_t r)
{
    m.newsize(n, k, r);
    my_func = NULL;
}

template <typename T>
Stencil3D<T>::Stencil3D(stencil3d_func f, integer_t n, integer_t k, integer_t r)
{
    m.newsize(n, k, r);
    my_func = f;
}

template <typename T>
void Stencil3D<T>::stencil_reg(stencil3d_func f)
{
    my_func = f;
}

template <typename T>
void Stencil3D<T>::stencil_reg(integer_t num, integer_t position[][3])
{
    integer_t move1 = m.dim2()*m.dim3();
    integer_t move2 = m.dim3();
    integer_t move3 = 1;
    Vector_hpc<T> sum(m.dim());


    for (integer_t i = 0; i<m.dim(); i++) {
        T tmp = 0.;
        for (integer_t k=0; k<num; k++) {
            tmp += *(m.ptr()+move1*position[k][0]+move2*position[k][1]+move3*position[k][2]);
        }
        sum[i] = tmp;
    }

    for (integer_t i = 0; i<m.dim(); i++) {
        *(m.ptr()+i) = (1.0/num)*sum[i];
    }
}


template <typename T>
void Stencil3D<T>::boundary(stencil3d_func f)
{
    if (m.dim() != 0) {
        f(m);
    }
}

template <typename T>
void Stencil3D<T>::run()
{
    if (my_func == NULL) {
        std::cerr << "Error: Stencil3D run has no stencil regedit! " << std::endl;
        return;
    }

    if (m.dim() != 0) {
        my_func(m);
    }
}


