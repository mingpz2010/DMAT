/*
    <Cmv.h: Basic abstract data type for DMAT project.>
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

#ifndef CMV_H_
#define CMV_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef integer_t
#undef integer_t
#endif

typedef long integer_t;

//
// Fortran changed to C/C++ basic data type
//
//      C++ type        name          Fortran type
// ----------------------------------------------
typedef unsigned char   BYTE;       // byte
typedef short int       INTEGER2;   // integer*2
typedef long            INTEGER;    // integer
typedef int             LOGICAL;    // logical
typedef float           REAL;       // real
typedef double          REAL8;      // real*8
typedef long double     REAL16;     // real*16

// support float, double, int and long basic data type (generic programming)
// " More C++ Idioms/Making New Friends "
// From the <More C++ Idioms> wiki book
// Declare friend function
template<typename T> class Vector_hpc;
template<typename T> Vector_hpc<T> operator+(const Vector_hpc<T> &c1, const Vector_hpc<T> &c2);
template<typename T> Vector_hpc<T> operator+(const Vector_hpc<T> &c1, T num);
template<typename T> Vector_hpc<T> operator+(T num, const Vector_hpc<T> &c1);
template<typename T> std::ostream& operator<<(std::ostream &s, const Vector_hpc<T> &A);
template<typename T> T max_NRM2(integer_t N, const Vector_hpc<T> &x1, const Vector_hpc<T> &x2);

template <typename T>
class Vector_hpc
{
protected:
    T *p_;
    integer_t dim_;
    int ref_;           // 0: own memory space; 1: point to another memory space
public:
    /*::::::::::::::::::::::::::*/
    /* Constructors/Destructors */
    /*::::::::::::::::::::::::::*/
    Vector_hpc() { p_ = NULL; dim_ = 0; ref_ = 0; }
    Vector_hpc(integer_t);
    Vector_hpc(integer_t, const T&);
    Vector_hpc(T*, integer_t);
    Vector_hpc(const T*, integer_t);
    Vector_hpc(const Vector_hpc &);
    ~Vector_hpc();
    /*::::::::::::::::::::::::::::::::*/
    /*  Indices and access operations */
    /*::::::::::::::::::::::::::::::::*/
    T& operator()(integer_t i) {
        return p_[i];
    }
    const T& operator()(integer_t i) const {
        return p_[i];
    }
    T& operator[](integer_t i) {
        return p_[i];
    }
    const T& operator[](integer_t i) const {
        return p_[i];
    }

    inline integer_t size() const { return dim_;}
    inline integer_t dim() const { return dim_;}
    inline integer_t ref() const { return ref_; }
    inline int null() const {return dim_== 0;}
    inline void zero() {
        for (integer_t i=0; i<dim_; i++) {
            p_[i] = 0;
        }
    }
    //
    // Create a new *uninitalized* vector of size N
    Vector_hpc<T> & newsize(integer_t);
    /*::::::::::::::*/
    /*  Assignment  */
    /*::::::::::::::*/
    Vector_hpc<T> & operator=(const Vector_hpc<T>&);
    Vector_hpc<T> & operator=(const T&);
    friend Vector_hpc<T> operator+ <>(const Vector_hpc<T> &c1, const Vector_hpc<T> &c2);
    friend Vector_hpc<T> operator+ <>(const Vector_hpc<T> &c1, T num);
    friend Vector_hpc<T> operator+ <>(T num, const Vector_hpc<T> &c1);

    // common functions
    void add(const Vector_hpc<T> &c1);
    void add(T *);
    void sub(const Vector_hpc<T> &c1);
    void sub(T *);
    void mul(T num);
    void div(T num);
    void fill(T num);
    T max();
    T min();
    T mean();

    // something related to Fortran
    void copyFortran(int ref, T *, INTEGER dim);

    friend std::ostream& operator<< <>(std::ostream &s, const Vector_hpc<T> &A);
    friend T max_NRM2 <>(integer_t N, const Vector_hpc<T> &x1, const Vector_hpc<T> &x2);
};

#include "Cmv.tpp"

#endif /* CMV_H_ */

