/*
	<vector_int.c, it is a vector ADT written by C language .>
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
#include "data.h"

static int 		*p_int 		= NULL;

int vector_int_alloc(unsigned int size, int *a)
{
	if (size <= 0) {
		return -1;
	}

	a = (int *)malloc(size * sizeof(int));

	return 0;
}

void vector_free(int *a)
{
	if (a == NULL) {
		return;
	}

	p_int = (int *)a;
	free(p_int);
	p_int = NULL;
}


