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

#define SIZE    10
#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test Vector_hpc class!\n");

    Vector_hpc<double> v1(SIZE);
    Vector_hpc<double> v2(SIZE);

    for (int i=0; i<SIZE; i++) {
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
    integer_t index;

    ans = v1.max();
    std::cout<< "v1 max : "<< ans << std::endl;
    ans = v1.mean();
    std::cout<< "v1 mean : "<< ans << std::endl;
    ans = v1.min();
    std::cout<< "v1 min : "<< ans << std::endl;

    // Basic arthmetic operation testing
    double tmp[SIZE] = {1, 2, 3, 4, 5, 6};
    v1.add(v2);
    v1.add(tmp);
    v1.sub(v2);
    v1.sub(tmp);
    v1.mul(3.14159);
    v1.div(2);
    v1.fill(5.0);
    std::cout << v1 << std::endl;
    std::cout << v1+v2 << std::endl;

    // BLAS operation testing
    v1.fill(2.0);
    v2.fill(3.0);
    v1.dswap(v1);
    v1.dswap(v2);
    v1.dscal(PI);
    v1.dcopy(v2);
    v1.daxpy(2.0, v2);
    ans = v1.ddot(v2);
    ans = v1.dnrm2();
    ans = v1.dasum();
    index = v1.idamax(ans);
    std::cout<<"v1 IDAMAX : index = "<<index<<" , Value = "<<ans<<std::endl;

    return 0;
}


