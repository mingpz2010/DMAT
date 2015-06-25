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

template <typename T>
Matrix_hpc<T>::Matrix_hpc() : Vector_hpc<T>()
{
	dim1_ = dim2_ = 0;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc(integer_t m, integer_t n)
{
	if (m<=0 || n<=0) {
	    Vector_hpc<T>();
	    dim1_ = dim2_ = 0;
	    std::cerr << "Error: bad value in Matrix_hpc constructor " << std::endl;
	    return;
	}
	dim1_ = m;
	dim2_ = n;
	Vector_hpc<T>(m*n);
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc<T> & newsize(integer_t m, integer_t n)
{
    if (m<0 || n<0) {
        return;
    }
    if (p_) delete [] p_;
    dim1_ = m;
    dim2_ = n;
    p_ = new T[m*n];
    if (p_ == NULL)
    {
    	std::cerr << "Error : NULL pointer in operator= Matrix_hpc newsize" << std::endl;
        exit(-1);
    }
    dim_ = m*n;

    return *this;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc<T> & operator=(const Matrix_hpc<T>& M)
{
    integer_t N = m.dim_;
    integer_t i;

    if (&M == this) {
        return *this;
    }

    // no need to test for overlap, since this region is new
    for (i =0; i< N; i++)       // careful not to use bcopy()
        p_[i] = m.p_[i];        // here, but double::operator= double.

    return *this;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc<T> & operator=(const T& m)
{
    // unroll loops to depth of length 8
    integer_t N = size();
    integer_t Nminus8 = N-8;
    integer_t i;

    for (i=0; i<Nminus8; )
    {
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
    }

    for (; i<N; p_[i++] = m);   // finish off last piece...

    return *this;
}

template <typename T>
void Matrix_hpc<T>::add(const Matrix_hpc<T> &c1)
{
    if (dim_ <= 0 || dim_ != c1.dim_) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] += c1.p_[i];
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const Matrix_hpc<T> &M)
{
    integer_t N = M.size();
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    for (integer_t i=0; i< i_max; i++) {
        for (integer_t j=0; j< j_max; j++) {
            s << M(i, j) << " ";
        }
    }
    
    s << std::endl;

    return s;
}


