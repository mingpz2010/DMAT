/*
	<vector_float.c, it is a vector ADT written by C language .>
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
#include <stdio.h>
#include "data.h"

static float 	*p_float 	= NULL;

int vector_float_alloc(unsigned int size, float *a)
{
	if (size <= 0) {
		return -1;
	}

	a = (float *)malloc(size * sizeof(float));

	return 0;
}

void vector_float_free(float *a)
{
	if (a == NULL) {
		return;
	}

	p_float = (float *)a;
	free(p_float);
	p_float = NULL;
}



