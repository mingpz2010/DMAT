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
#include <cstdlib>
#include <omp.h>

#define BLOCKSIZE   100
#define N           1000

using namespace std;

void func1()
{
    int nthreads, tid;

    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        printf("HELLO WORLD from thread = %d\n", tid);
        if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("NUMBER of threads = %d\n", nthreads);
        }
    }
}

void func2()
{
    int i, size;
    double a[N], b[N], c[N];

    for (int i=0; i<N; i++) {
        a[i] = b[i] = i * 1.0;
    }
    size = BLOCKSIZE;

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(dynamic,size) nowait
        for(int i=0; i<N; i++) {
            c[i] = a[i] + b[i];
        } /* End of parallel section */
    }

    printf("%.6lf + %.6lf = %.6lf\n", a[N/2], b[N/2], c[N/2]);
}

void func3()
{
    double a[N], b[N], c[N], d[N];

    for (int i=0; i<N; i++) {
        a[i] = i * 1.5;
        b[i] = i + 22.35;
    }

    #pragma omp parallel default(shared)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            for (int i=0; i<N; i++) {
                c[i] = a[i] + b[i];
            }
            #pragma omp section
            for (int i=0; i<N; i++) {
                d[i] = a[i] * b[i];
            }
        }  /* End of sections */
    }  /* End of parallel section */

    printf("%.6lf + %.6lf = %.6lf\n", a[N/2], b[N/2], c[N/2]);
    printf("%.6lf * %.6lf = %.6lf\n", a[N/2], b[N/2], d[N/2]);
}

void func4()
{
    int x;

    x = 0;
    #pragma omp parallel shared(x)
    {
        #pragma omp critical
        x += 1;
    }  /* End of parallel section */
}

void func5()
{
    int i, n, size;
    double a[100], b[100], result;

    n = 100;
    size = 10;
    result = 0.;
    for (int i=0; i<n; i++) {
        a[i] = i * 1.0;
        b[i] = i * 2.0;
    }

    for (i = 0; i<n; i++) {
        result += a[i] * b[i];
    }
    printf("RESULT SERIAL = %.6lf\n", result);

    result = 0.;
    #pragma omp parallel for \
		default(shared) private(i) \
		schedule(static,size) \
		reduction(+:result)
    for (i = 0; i<n; i++) {
        result += a[i] * b[i];
    }

    printf("FINAL RESULT = %.6lf\n", result);
}

int main(int argc, char *argv[])
{
    func1();
    func2();
    func3();
    func4();
    func5();

    return 0;
}


