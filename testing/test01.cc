/*
	<test01.cc, testing for c++ language library version.>
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
#include "../src/Cmv.h"

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test Vector_hpc class and Matrix_hpc class!\n");

    Vector_hpc<double> v1(10);
    Vector_hpc<double> v2(10);

    for (int i=0; i<10; i++) {
        v1(i) = i*PI;
    }
    std::cout << "output v1:"<< std::endl;
    std::cout << v1 << std::endl;
    std::cout << "output v2:"<< std::endl;
    std::cout << v2 << std::endl;
    std::cout << "output v2:"<< std::endl;
    v2 = v1;
    std::cout << v2 << std::endl;

    double ans;

    ans = v1.max();
    std::cout<< "v1 max : "<< ans << std::endl;
    ans = v1.mean();
    std::cout<< "v1 mean : "<< ans << std::endl;
    ans = v1.min();
    std::cout<< "v1 min : "<< ans << std::endl;
    v1.mul(3.14159);
    std::cout << v1 << std::endl;

    std::cout << v1+v2 << std::endl;

    return 0;
}


