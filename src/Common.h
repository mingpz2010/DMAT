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
#ifndef COMMON_H_
#define COMMON_H_

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

typedef enum Matrix_store_manner {
    ROW_MAJOR,
    COL_MAJOR
}matrix_manner_t;

typedef enum Matrix_property {
    M_NORMAL,
    M_BANDED,
    M_SYMMETRIC,
    M_SYMMETRIC_BANDED,
    M_TRIANGULAR
}matrix_property_t;


#endif /* COMMON_H_ */
