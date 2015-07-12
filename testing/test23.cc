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
#include "../src/Dimscal.h"

#define SIZEA   100
#define SIZEB   200
#define SIZEC   2000

double now()
{
    std::clock_t t = clock();

    return static_cast<double>(t)/CLOCKS_PER_SEC;
}

int main(int argc, char *argv[])
{
    double start, end;
    double a = 0.9, b = 1e-10, c = 2.5;
    Dimscal<double> x(SIZEA, SIZEB, SIZEC);

    for (int i=0; i<SIZEA; i++) {
        for (int j=0; j<SIZEB; j++) {
            for (int k=0; k<SIZEC; k++) {
                x(i,j,k) = random();
            }
        }
    }

    start = now();
    x.blas_op(a,b,c);
    end = now();

    // std::cout<<"Operation 2(arithmetic) cost time "<< end-start <<" (s)" << std::endl;
    printf("Operation 2(arithmetic) cost time %.8lf\n", end-start);

    return 0;
}


