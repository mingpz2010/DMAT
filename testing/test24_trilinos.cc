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

#include <iostream>
#include <cstdio>
#include <ctime>
#include "../src/Stencil3D.h"

#define SIZEA   10
#define SIZEB   20
#define SIZEC   5

double now()
{
    std::clock_t t = clock();

    return static_cast<double>(t)/CLOCKS_PER_SEC;
}

template <typename T>
void stencil_calc(Dimscal<T>& d)
{
    // Original Manner
    for(int i=0; i<d.dim1(); i++) {
        for (int j=0; j<d.dim2(); j++) {
            for (int k=0; k<d.dim3(); k++) {
                d(i,j,k) = (1.0/6)*(d(i-1,j,k)+d(i+1,j,k)+d(i,j+1,k)+
                        d(i,j-1,k)+d(i,j,k+1)+d(i,j,k-1));
            }
        }
    }
}

int main(int argc, char *argv[])
{
    double start, end;
    Stencil3D stencil<double>(SIZEA, SIZEB, SIZEC);
    long pos[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};

#if 0
    // Original Manner
    stencil.stencil_reg(stencil_calc);
#endif

    // New Manner
    stencil.stencil_reg(6, pos);

    start = now();
    stencil.run();
    end = now();

    std::cout<<"Operation 3(Stencil3D) cost time "<< end-start <<" (s)" << std::endl;

    return 0;
}



