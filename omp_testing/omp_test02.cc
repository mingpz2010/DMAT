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
#include <complex>
#include <cmath>
#include <omp.h>

using namespace std;

typedef std::complex<double> my_complex;

void demo1()
{
    const int size = 256;
    double sinTable[size];
    printf("Demo1:\n");
    for (int n=0; n<size; ++n) {
        sinTable[n] = std::sin(2 * M_PI * n / size);
    }
}

int MandelbrotCalculate(my_complex c, int maxiter)
{
    my_complex z = c;
    int n = 0;
    for (; n<maxiter; ++n) {
        if (std::abs(z) >= 2.0)  break;
        z = z * z + c;
    }

    return n;
}

void demo2()
{
    printf("Demo2:\n");
    const int width = 78, height = 44, num_pixels = width*height;
    const my_complex center(-.7, 0), span(2.7, -(4/3.0)*2.7*height/width);
    const my_complex begin = center-span/2.0, end = center+span/2.0;
    const int maxiter = 100000;

    #pragma omp parallel for ordered schedule(dynamic)
    for(int pix=0; pix<num_pixels; ++pix)
    {
        const int x = pix%width, y = pix/width;

        my_complex c = begin + my_complex(x * span.real() / (width +1.0),
                y * span.imag() / (height+1.0));

        int n = MandelbrotCalculate(c, maxiter);
        if(n == maxiter) n = 0;
        #pragma omp ordered
        {
            char c = ' ';
            if(n > 0)
            {
                static const char charset[] = ".,c8M@jawrpogOQEPGJ";
                c = charset[n % (sizeof(charset)-1)];
            }
            std::putchar(c);
            if(x+1 == width) std::puts("|");
        }
    }

    #pragma omp barrier
}

void report_num_threads(int level)
{
    #pragma omp single
    {
        printf("Level %d: number of threads in the team - %d\n",
                level, omp_get_num_threads());
    }
}

void demo3()
{
    printf("Demo3:\n");

    omp_set_dynamic(0);
#if 0
    printf("NESTED FLAG = %d\n", omp_get_nested());
    omp_set_nested(1);
    printf("NEW NESTED FLAG = %d\n", omp_get_nested());
#endif

    printf("SB\n");

    #pragma omp parallel num_threads(4)
    {
        report_num_threads(1);

        #pragma omp parallel num_threads(4)
        {
            report_num_threads(2);
        }
    }  /* End of PARALLEL section */
}

void demo4()
{
    printf("Demo4:\n");
    omp_set_dynamic(0);
    printf("OUTPUT From FIRST LOOP\n");
    #pragma omp parallel num_threads(4)
    {
        report_num_threads(1);
        #pragma omp parallel num_threads(4)
        {
            report_num_threads(2);
        }
    }
    printf("OUTPUT From SECOND LOOP\n");
    #pragma omp parallel num_threads(8)
    report_num_threads(1);
}

void report_num_threads_ForLoop(int id, int level)
{
    #pragma omp single
    {
        printf("Level %d[%d]: number of threads in the team - %d\n",
                level, id, omp_get_num_threads());
    }
}

void demo5()
{
    printf("Demo5:\n");
    omp_set_dynamic(0);
    #pragma omp parallel for num_threads(2)
    for (int i=101; i<105; i++) {
        report_num_threads_ForLoop(i, 1);
        #pragma omp parallel num_threads(3)
        {
            report_num_threads_ForLoop(i, 2);
        }
    }
}

void demo6()
{
    double a[4], b[4];

    printf("Demo6:\n");
    for (int i=0; i<4; i++) {
        a[i] = b[i] = i + 2;
    }
    #pragma omp parallel for ordered schedule(dynamic)
    for(int i=0; i < 4; ++i) {
        a[i] *= (b[i]/3);
        #pragma omp ordered
        printf("i = %d, a[%d] = %.6lf\n", i, i, a[i]);
    }
}

int main(int argc, char *argv[])
{
    printf("\nEXPLAIN THE OPENMP in C++\n\n");
    demo1();
    // demo2();
    demo6();
    demo3();
    demo4();
    demo5();

    return 0;
}

