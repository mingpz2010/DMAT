/*
 * Cmv.cc
 *
 *  Created on: 2015.6.11
 *      Author: PingzhouMing
 */
#include <iostream>
#include <fstream>
#include "Cmv.h"

using namespace std;

Vector_double::Vector_double() : p_(0), dim_(0) {};

Vector_double::Vector_double(integer_t n) : p_(new double[n]), dim_(n)
{
    if (p_ == NULL) {
        cerr << "Error: NULL pointer in Vector_double(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(-1);
    }
}

Vector_double::Vector_double(integer_t n, const double& v) :
        p_(new double[n]), dim_(n)
{
    if (p_ == NULL)
    {
        cerr << "Error: NULL pointer in Marray_1D_double(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(-1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = v;
}

Vector_double::Vector_double(double* d, integer_t n) : p_(new double[n]),
      dim_(n)
{
    if (p_ == NULL)
    {
        cerr << "Error: Null pointer in Marray_1D_double(double*, int) " << endl;
        exit(1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = d[i];

}

Vector_double::Vector_double(const double* d, integer_t n) : p_(new double[n]),
      dim_(n)
{
    if (p_ == NULL)
    {
        cerr << "Error: Null pointer in Marray_1D_double(double*, int) " << endl;
        exit(-1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = d[i];

}

Vector_double::Vector_double(const Vector_double & m) : p_(new double[m.dim_]),
    dim_(m.dim_)
{
    if (p_ == NULL)
    {
        cerr << "Error:  Null pointer in Marray_1D_double(const Marray_1D_double&); " << endl;
        exit(-1);
    }

    integer_t N = m.dim_;

    for (integer_t i=0; i<N; i++)
        p_[i] = m.p_[i];
}

Vector_double::~Vector_double()
{
    if (p_) delete [] p_;
}


Vector_double& Vector_double::newsize(integer_t n)
{
    if (dim_ != n )                     // only delete and new if
    {                                   // the size of memory is really
        if (p_) delete [] p_;           // changing, otherwise just
        p_ = new double[n];              // copy in place.
        if (p_ == NULL)
        {
            cerr << "Error : NULL pointer in operator= newsize" << endl;
            exit(-1);
        }
        dim_ = n;
    }

    return *this;
}

Vector_double& Vector_double::operator=(const double & m)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
    }

    for (; i<N; p_[i++] = m);   // finish off last piece...

    return *this;
}

Vector_double& Vector_double::operator=(const Vector_double & m)
{

    integer_t N = m.dim_;
    integer_t i;

    newsize(N);

    // no need to test for overlap, since this region is new
    for (i =0; i< N; i++)       // careful not to use bcopy()
        p_[i] = m.p_[i];        // here, but double::operator= double.

    return *this;
}

Vector_double operator+(const Vector_double &c1, const Vector_double &c2)
{
    Vector_double c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] = c1[i] + c2[i];
    }

    return c;
}

Vector_double operator+(const Vector_double &c1, double num)
{
    Vector_double c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] += num;
    }

    return c;
}

Vector_double operator+(double num, const Vector_double &c1)
{
    Vector_double c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] += num;
    }

    return c;
}

void Vector_double::add(const Vector_double &c1)
{
    if (dim_ <= 0 || dim_ != c1.dim_) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] += c1[i];
    }
}

void Vector_double::add(double *m)
{
    if (dim_ <= 0) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] += m[i];
    }
}

void Vector_double::sub(const Vector_double &c1)
{
    if (dim_ <= 0 || dim_ != c1.dim_) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] -= c1[i];
    }
}

void Vector_double::sub(double *m)
{
    if (dim_ <= 0) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] -= m[i];
    }
}

void Vector_double::mul(double num)
{
    if (dim_ <= 0 )  return;

    integer_t i;

    for (i=0; i<dim; i++) {
        p_[i] *= num;
    }
}

void Vector_double::div(double num)
{
    if (dim_ <= 0 )  return;

    if (num == 0) {
        cerr << "div 0 error, so we don't process it"<< endl;
        return;
    }

    integer_t i;

    for (i=0; i<dim; i++) {
        p_[i] /= num;
    }
}

double Vector_double::max()
{
    if (dim_ <= 0) { return 0.; }
    if (dim_ == 1) { return p_[0]; }

    double ans;
    integer_t i;
    ans = p_[0];
    for (i=1; i<dim; i++) {
        if( ans < p_[i]) ans = p_[i];
    }
    return ans;
}

double Vector_double::min()
{
    if (dim_ <= 0) { return 0.; }
    if (dim_ == 1) { return p_[0]; }

    double ans;
    integer_t i;
    ans = p_[0];
    for (i=1; i<dim; i++) {
        if( ans > p_[i]) ans = p_[i];
    }
    return ans;
}

double Vector_double::mean()
{
    if (dim_ <= 0)  return 0.;

    double ans;
    integer_t i;
    for ( ans = 0., i=0; i<dim; i++) {
        ans += p_[i];
    }

    return ans/dim_;
}

void Vector_double::copyFortran(int ref, REAL8 *from, INTEGER dim)
{
    ref_ = ref;

    if (p_ != NULL) {
        delete []p_;
    }
    dim_ = dim;
    if (ref_ == 0) {
        p_ = new double[dim_];
    } else {
        p_ = from;
    }
}

ostream& operator<<(ostream& s, const Vector_double& V)
{
    integer_t N = V.size();

    for (integer_t i=0; i< N; i++)
        s << V(i) << endl;

    return s;
}

Face_current::Face_current(integer_t dim1, integer_t dim2, integer_t dim3, integer_t dim4)
                    : Vector_double(dim1*dim2*dim3*dim4)
{
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    dim_ = dim1*dim2*dim3*dim4;
    w1 = dim2_*dim3_*dim4_;
    w2 = dim3_*dim4_;
    w3 = dim4_;
}

Dimensional_scal::Dimensional_scal(integer_t dim1, integer_t dim2, integer_t dim3)
                    : Vector_double(dim1*dim2*dim3)
{
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim_ = dim1*dim2*dim3;
    w1 = dim2*dim3;
    w2 = dim3;

    p_ = new double[dim_];
}

double Dimensional_scal::maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3)
{
    if (dim_ <= 0) {
        dim1 = dim2 = dim3 = -1;
        return 0.;
    }
    if (dim_ == 1) {
        dim1 = dim2 = dim3 = 0;
        return p_[0];
    }

    integer_t i,j,k;
    integer_t index;
    double ans;

    ans = 0.;
    for (i=0; i<dim1_; i++) {
        index = i*w1;
        for (j=0; j<dim2_; j++) {
            index += j*w2;
            for (k=0; k<dim3_; k++) {
                index += k;
                if (ans < p_[index] ) {
                    ans = p_[index];
                    dim1 = i; dim2 = j; dim3 = k;
                }
            }
        }
    }

    return ans;
}

void Dimensional_scal::copyFortran(int ref, REAL8 *from, INTEGER dim1,
        INTEGER dim2, INTEGER dim3)
{
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim_ = dim1*dim2*dim3;
    w1 = dim2*dim3;
    w2 = dim3;

    Vector_double::copyFortran(ref, from, dim1*dim2*dim3);
}

