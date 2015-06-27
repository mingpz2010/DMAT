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
Flux<T>::Flux(integer_t m, integer_t n, integer_t k)
{
	if (m<=0 || n<=0 || k<=0) {
	    Vector_hpc<T>::p_ = NULL;
	    Vector_hpc<T>::dim_ = 0;
	    dim1_ = dim2_ = dim3_ = 0;
	    std::cerr << "Error: bad value in Flux constructor " << std::endl;
	    return;
	}
	dim1_ = m;
	dim2_ = n;
	dim3_ = k;
	w1 = n*k;
	w2 = k;
	Vector_hpc<T>::p_ = new T[m*n*k];
	Vector_hpc<T>::dim_ = m*n*k;
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Flux<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
}

template <typename T>
Flux<T>& Flux<T>::newsize(integer_t m, integer_t n, integer_t k)
{
    if (m<=0 || n<=0 || k<=0) {
        return;
    }
    if (Vector_hpc<T>::p_) delete [] Vector_hpc<T>::p_;
    dim1_ = m;
    dim2_ = n;
    dim3_ = k;
    w1 = n*k;	w2 = k;
    Vector_hpc<T>::p_ = new T[m*n*k];
    if (Vector_hpc<T>::p_ == NULL)
    {
    	std::cerr << "Error : NULL pointer in operator= Flux newsize" << std::endl;
        exit(-1);
    }
    Vector_hpc<T>::dim_ = m*n*k;

    return *this;
}

template <typename T>
Flux<T>& Flux<T>::operator=(const Flux<T>& M)
{
    integer_t N = M.size();
    integer_t i;

    if (&M == this) {
        return *this;
    }

    // no need to test for overlap, since this region is new
    for (i =0; i< N; i++)       // careful not to use bcopy()
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];        // here, but double::operator= double.

    return *this;
}

template <typename T>
Flux<T>& Flux<T>::operator=(const T& m)
{
    // unroll loops to depth of length 8
    integer_t N = size();
    integer_t Nminus8 = N-8;
    integer_t i;

    for (i=0; i<Nminus8; )
    {
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
        Vector_hpc<T>::p_[i++] = m;
    }

    for (; i<N; Vector_hpc<T>::p_[i++] = m);   // finish off last piece...

    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const Flux<T> &M)
{
    integer_t N = M.size();
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();
    integer_t k_max = M.dim3();

    for (integer_t i=0; i< i_max; i++) {
        for (integer_t j=0; j< j_max; j++) {
        	for (integer_t k=0; k< k_max; k++) {
            	s << M(i, j, k) << std::endl;
            }
        }
    }
    
    s << std::endl;

    return s;
}

