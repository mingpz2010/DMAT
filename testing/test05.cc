/*
    <test05.cc, testing for c++ language library version.>
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
#include "../src/Flux.h"

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test Flux class!\n");

    Flux<double> m1(2, 2, 2, 10);
    Flux<double> m2(2, 2, 2, 10);

    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
            for (int k=0; k<2; k++) {
                for (int r=0; r<10; r++) {
                    m1(i, j, k, r) = (i+j+k+r)*PI;
                }
            }
        }
    }

    std::cout << "output m1:"<< std::endl;
    std::cout << m1 << std::endl;
    std::cout << "output m2:"<< std::endl;
    std::cout << m2 << std::endl;
    std::cout << "output m2:"<< std::endl;
    m2 = m1;
    std::cout << m2 << std::endl;

    m1 = 1e-14;
    std::cout << m1 << std::endl;

    return 0;
}



