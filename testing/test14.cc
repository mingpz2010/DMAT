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
    TRACE_PRINT("start to test SparseMatrix TDS manner rightness?\n");

    Matrix_hpc<double> m_diagonal(6, 6);
    m_diagonal(0,0) = 10; m_diagonal(0,1) = -2;
    m_diagonal(1,1) = 9; m_diagonal(1,0) = 3; m_diagonal(1,2) = 6;
    m_diagonal(2,2) = 8; m_diagonal(2,1) = 7; m_diagonal(2,3) = 7;
    m_diagonal(3,3) = 7; m_diagonal(3,2) = 8; m_diagonal(3,4) = 5;
    m_diagonal(4,4) = 9; m_diagonal(4,3) = 9; m_diagonal(4,5) = 13;
    m_diagonal(5,5) = -1; m_diagonal(5,4) = 2;
    SparseMatrix<double> m1(TDS_MANNER, m_diagonal);

    std::cout << "output m1:"<< std::endl;
    std::cout << m1 << std::endl;

    double x[6];
    double b[6] = {1, 1, 1, 1, 1, 1};

    m1.chase_method(6, x, b);
    std::cout << "output chase method solution:" << std::endl;
    for (int i=0; i<6; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "standard solution are:" << std::endl;
    std::cout << "0.04516851 -0.27415747 0.23599264 0.22018732 -0.48268941 0.25865529" << std::endl;

    TRACE_PRINT("start to test SparseMatrix class efficiency!\n");

    return 0;
}

