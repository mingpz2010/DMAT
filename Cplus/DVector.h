/*
	<DVector.h, it is a vector template which is used for reactor calculation.>
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

#ifndef DVECTOR_H_
#define DVECTOR_H_

#include <iostream>

template <typename T>
class DVector
{
    protected:
        T *p_;
        unsigned int dim_;
    public:
        DVector() { p_ = NULL; dim_ = 0; }
        DVector(unsigned int n) {
            if (n <= 0) { p_ = NULL; dim_ = 0; }
            p_ = new T[n]; dim_ = n;
        }

        T& operator()(unsigned int i) {
            return p_[i];
        }
        const T& operator()(unsigned int i) const {
            return p_[i];
        }
        T& operator[](unsigned int i) {
            return p_[i];
        }
        const T& operator[](unsigned int i) const {
            return p_[i];
        }

        inline unsigned int size() const { return dim_; }
        inline unsigned int dim() const { return dim_; }
        inline int null() const { return dim_== 0; }

        DVector& operator=(const DVector&);
        DVector& operator=(const T&);
};


#endif /* DVECTOR_H_ */

