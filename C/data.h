/*
	<data.h, it is a header file .>
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

#ifndef DATA_H_
#define DATA_H_

typedef enum data_type {
    CHAR,
    INT,
    LONG,
    FLOAT,
    DOUBLE
}data_type_t;

int vector_double_alloc(unsigned int size, double *a);
void vector_double_free(double *a);

#endif /* DATA_H_ */

