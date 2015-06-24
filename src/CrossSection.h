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

#ifndef CROSSSECTION_H_
#define CROSSSECTION_H_

#include "Cmv.h"

// It derived from Vector_hpc, 2D ADT with { energy group + value }
template<typename T>
class CrossSection : public Vector_hpc<T>
{

};

#include "CrossSection.tpp"

#endif /* CROSSSECTION_H_ */