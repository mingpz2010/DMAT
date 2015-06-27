/*
    <test06.cc, testing for c++ language library version.>
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
#include "../src/Flux.h"

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

double now()
{
    std::clock_t t = clock();

    return static_cast<double>(t)/CLOCKS_PER_SEC;
}

void Flux_benchmark(long size)
{
    Flux<double> v1(size, size, size);
    Flux<double> v2(size, size, size);

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            for (int k=0; k<size; k++) {
                v1(i,j,k) = (i+j+k)*PI;
            }
        }
    }
    v2 = v1;
}

void basic_benchmark(long size)
{
    double*** v1 = new double**[size];
    for(int i = 0; i < size; ++i) {
        v1[i] = new double*[size];
        for (int j=0; j<size; ++j) {
            v1[i][j] = new double[size];
        }
    }
    double*** v2 = new double**[size];
    for(int i = 0; i < size; ++i) {
        v2[i] = new double*[size];
        for (int j=0; j<size; ++j) {
            v2[i][j] = new double[size];
        }
    }

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            for (int k=0; k<size; k++) {
                v1[i][j][k] = (i+j+k)*PI;
            }
        }
    }

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            for (int k=0; k<size; k++) {
                v2[i][j][k] = v1[i][j][k];
            }
        }
    }

    for(int i = 0; i < size; ++i) {
        for (int j=0; j<size; ++j) {
            delete[] v2[i][j];
            delete[] v1[i][j];
        }
        delete[] v2[i];
        delete[] v1[i];
    }
}

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to test Matrix_hpc class efficiency!\n");

    double start, end;

    start = now();
    Matrix_hpc_benchmark(1000);
    end = now();

    TRACE_PRINT("[%d] Matrix_hpc<double> run time is %.6lf sec.\n", 1000, end-start);

    start = now();
    basic_benchmark(1000);
    end = now();

    TRACE_PRINT("[%d] basic array[][] run time is %.6lf sec.\n", 1000, end-start);

    return 0;
}


