/*
    <test03.cc, testing for c++ language library version.>
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
#include "../src/Matrix_hpc.h"

#define SIZE    8
#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test Matrix_hpc class!\n");

    Matrix_hpc<double> m1(SIZE, SIZE);
    Matrix_hpc<double> m2(SIZE, SIZE);

    for (int i=0; i<SIZE; i++) {
        for (int j=0; j<SIZE; j++) {
            m1(i, j) = (i+j)*PI;
        }
    }

    std::cout << "output m1:"<< std::endl;
    std::cout << m1 << std::endl;
    std::cout << "output m2:"<< std::endl;
    std::cout << m2 << std::endl;
    std::cout << "output m2:"<< std::endl;
    m2 = m1;
    std::cout << m2 << std::endl;
    m1.add(m2);
    std::cout << m1 << std::endl;
    m1.newsize(1, 1);
    std::cout << m1 << std::endl;

    // BLAS operation testing
    double ans;
    Vector_hpc<double> v1(size);
    Vector_hpc<double> v2(size);
    for (int i=0; i<SIZE; i++) {
        v1[i] = (i+1) * 2.5;
        v2[i] = 0.;
    }

    m1.dswap(m2);
    m1.dscal(PI);
    m1.dcopy(m2);
    m1.daxpy(0.2, m2);
    ans = m1.dnrm2();
    ans = m1.dasum();
    m1.dgemv(v1);
    m1.dgemv(v1,v2);

    return 0;
}





