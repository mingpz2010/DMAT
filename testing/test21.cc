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
#include <cmath>
#include "../src/Cmv.h"
#include "../src/SparseMatrix.h"

typedef struct material {
    double tr[2];
    double a[2];
    double vf[2];
    double kf[2];
    double tr1to2;
    double d1;
    double d2;
    double exposure;
}material_t;

material_t mat1;
material_t mat2;
material_t mat3;

using namespace std;

void diffusion_solver(int single_mesh)
{
    int geom[13] = {3,3,3,3,3,3,3,3,3,3,3,3,3};
    double chi1, chi2, mesh_space;
    long int mesh;
    Vector_hpc<double> flux1;
    Vector_hpc<double> flux2;
    Vector_hpc<double> source1;
    Vector_hpc<double> source2;
    Vector_hpc<double> p;
    Vector_hpc<double> flux1_last;
    Vector_hpc<double> flux2_last;
    SparseMatrix<double> mat_of_coeff1;
    SparseMatrix<double> mat_of_coeff2;
    double keff, keff_last, cjk, cjf1, cjf2;
    material_t *mat;

    mat1.exposure = 0.2;
    mat1.tr[0]=2.368355e-1; mat1.tr[1]=9.082422e-1;
    mat1.a[0]=8.603111e-3; mat1.a[1]=7.853449e-2;
    mat1.vf[0]=6.160544e-3; mat1.vf[1]=1.207603e-1;
    mat1.kf[0]=8.158099e-14; mat1.kf[1]=1.599168e-12;
    mat1.tr1to2=1.708253e-2;
    mat1.d1 = 1./(3*mat1.tr[0]);
    mat1.d2 = 1./(3*mat1.tr[1]);
    cout << "MAT1 : " << mat1.d1 <<" , " << mat1.d2 << endl;

    mat2.exposure = 8.11;
    mat2.tr[0]=2.367121e-1; mat2.tr[1]=9.239672e-1;
    mat2.a[0]=9.150915e-3; mat2.a[1]=8.516051e-2;
    mat2.vf[0]=5.585696e-3; mat2.vf[1]=1.250261e-1;
    mat2.kf[0]=7.144699e-14; mat2.kf[1]=1.599212e-12;
    mat2.tr1to2=1.690581e-2;
    mat2.d1 = 1./(3*mat2.tr[0]);
    mat2.d2 = 1./(3*mat2.tr[1]);
    cout << "MAT2 : " << mat2.d1 <<" , " << mat2.d2 << endl;

    mat3.exposure = 16.55;
    mat3.tr[0]=2.366212e-1; mat3.tr[1]=9.308326e-1;
    mat3.a[0]=9.668583e-3; mat3.a[1]=8.506164e-2;
    mat3.vf[0]=5.050670e-3; mat3.vf[1]=1.188626e-1;
    mat3.kf[0]=6.310672e-14; mat3.kf[1]=1.485153e-12;
    mat3.tr1to2=1.675986e-2;
    mat3.d1 = 1./(3*mat3.tr[0]);
    mat3.d2 = 1./(3*mat3.tr[1]);
    // cout << "MAT3 : " << mat3.d1 <<" , " << mat3.d2 << endl;
    printf("MAT3 : %.12lf , %.12lf\n", mat3.d1, mat3.d2);

    mat = &mat3;
    chi1 = 1.0; chi2 = 0.;
    keff = 1.0;
    mesh = single_mesh * 13; mesh_space = 20.0/single_mesh;
    flux1.newsize(mesh);
    flux2.newsize(mesh);
    source1.newsize(mesh);
    source2.newsize(mesh);
    p.newsize(mesh);
    flux1_last.newsize(mesh);
    flux2_last.newsize(mesh);

    cout << "1D Reactor geometry is : " << mesh << endl;

    double *dia = new double[mesh];
    double *left = new double[mesh];
    double *right = new double[mesh];

    int k;
    for (int i=0; i<13; i++) {
        k = 0;
        int idx = i * single_mesh;
        while (k < single_mesh) {
            dia[idx+k] = (2.*mat->d1)/(mesh_space*mesh_space)+mat->a[0]+mat->tr1to2;
            k = k + 1;
        }
    }
    for (int i=0; i<mesh; i++) {
        if (i>0) {
            left[i] = -mat->d1/(mesh_space*mesh_space);
        }
        if (i<mesh-1) {
            right[i] = -mat->d1/(mesh_space*mesh_space);
        }
    }
    left[0] = 0; right[mesh-1] = 0;
    dia[mesh-1] = 1;
    mat_of_coeff1.tds_alloc(mesh, dia, left, right);

#if 0
    cout << "mat_of_coeff1(DIA) : " << endl;
    printf("%.16lf %.16lf %.16lf", dia[0], dia[mesh/2], dia[mesh-1]);
    printf("\n");
    cout << "mat_of_coeff1(left) : " << endl;
    printf("%.16lf %.16lf %.16lf", left[0], left[mesh/2], left[mesh-1]);
    printf("\n");
    cout << "mat_of_coeff1(right) : " << endl;
    printf("%.16lf %.16lf %.16lf", right[0], right[mesh/2], right[mesh-1]);
    printf("\n");
#endif

    for (int i=0; i<13; i++) {
        k = 0;
        int idx = i * single_mesh;
        while (k < single_mesh) {
            dia[idx+k] = (2.*mat->d2)/(mesh_space*mesh_space)+mat->a[1];
            k = k + 1;
        }
    }
    for (int i=0; i<mesh; i++) {
        if (i>0) {
            left[i] = -mat->d2/(mesh_space*mesh_space);
        }
        if (i<mesh-1) {
            right[i] = -mat->d2/(mesh_space*mesh_space);
        }
    }
    left[0] = 0; right[mesh-1] = 0;
    dia[mesh-1] = 1;
    mat_of_coeff2.tds_alloc(mesh, dia, left, right);

    keff = 0.9;
    int itr = 1;
    for (int i=0; i<mesh; i++) {
        flux1(i) = keff/((mat->vf[0]+mat->vf[1])*(mesh*mesh_space));
        flux2(i) = flux1(i);
        p(i) = mat->vf[0]*flux1(i)+mat->vf[1]*flux2(i);
    }

    cjk = 1;
    cjf1 = 1;
    cjf2 = 1;
    cout << "Start to Iteration!" << endl;
    while (cjk >= 1e-16 || cjf1>=1e-16 || cjf2>=1e-16) {
        for (int i=0; i<mesh; i++) {
            source1(i) = (chi1*p(i))/keff;
        }
        flux1_last = flux1;
        flux2_last = flux2;
        source1(mesh-1) = 0.;
#if 0
        cout << "source1 : " << endl;
        printf("%.16lf %.16lf %.16lf", source1[0], source1[mesh/2], source1[mesh-1]);
        printf("\n");
#endif
        mat_of_coeff1.chase_method(mesh, flux1, source1);
#if 1
        cout << "flux1 : " << endl;
        printf("%.16lf %.16lf %.16lf", flux1[0], flux1[mesh/2], flux1[mesh-1]);
        printf("\n");
#endif
        for (int i=0; i<mesh; i++) {
            source2(i) = (chi2*p(i))/keff + mat->tr1to2 * flux1(i);
        }
        mat_of_coeff2.chase_method(mesh, flux2, source2);

        keff_last = keff;
        keff = 0;
        for (int i=0; i<mesh; i++) {
            p(i) = mat->vf[0]*flux1(i)+mat->vf[1]*flux2(i);
            keff += p(i) * mesh_space;
        }
        cjk = fabs(keff-keff_last)/keff_last;
        itr = itr + 1;
        cjf1 = max_NRM2<double>(mesh, flux1, flux1_last);
        cjf2 = max_NRM2<double>(mesh, flux2, flux2_last);
        if (itr == 2) {
            exit(0);
        }
        printf("itr = %d: keff = %.8lf, CJK = %.4le, CJF1 = %.4le, CJF2 = %.4le\n",
                itr-1, keff, cjk, cjf1, cjf2);
    }

    cout << "keff = " << keff << endl;

    if (dia != NULL) {
        delete[] dia;
    }
    if (left != NULL) {
        delete[] left;
    }
    if (right != NULL) {
        delete[] right;
    }
}

int main(int argc, char *argv[])
{
    std::cout << "Start to perform 1D diffusion calculation"<< std::endl;
    std::cout << "-----------------------------------------"<< std::endl;

    diffusion_solver(2);

    return 0;
}

