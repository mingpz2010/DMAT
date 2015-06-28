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
#include "../src/SparseMatrix.h"

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test SparseMatrix class!\n");

    Matrix_hpc<double> m1(6, 6);
    m1(0,0) = 10; m1(0,4) = -2;
    m1(1,0) = 3; m1(1,1) = 9; m1(1,5) = 3;
    m1(2,1) = 7; m1(2,2) = 8; m1(2,3) = 7;
    m1(3,0) = 3; m1(3,2) = 8; m1(3,3) = 7; m1(3,4) = 5;
    m1(4,1) = 8; m1(4,3) = 9; m1(4,4) = 9; m1(4,5) = 13;
    m1(5,1) = 4; m1(5,4) = 2; m1(5,5) = -1;

    SparseMatrix<double> m2(m1);
    SparseMatrix<double> m3(CCS_MANNER, m1);

    std::cout << "output m1:"<< std::endl;
    std::cout << m1 << std::endl;
    std::cout << "output m2:"<< std::endl;
    std::cout << m2 << std::endl;
    std::cout << "output m3:"<< std::endl;
    std::cout << m3 << std::endl;

    return 0;
}

