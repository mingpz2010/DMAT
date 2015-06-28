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
Matrix_hpc<T>::Matrix_hpc()
{
    type = 0;
    Vector_hpc<T>::p_ = NULL;
    Vector_hpc<T>::dim_ = 0;
	dim1_ = dim2_ = 0;
	nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc(integer_t m, integer_t n)
{
	if (m<=0 || n<=0) {
	    Vector_hpc<T>::p_ = NULL;
	    Vector_hpc<T>::dim_ = 0;
	    dim1_ = dim2_ = 0;
	    std::cerr << "Error: bad value in Matrix_hpc constructor " << std::endl;
	    return;
	}
	type = 0;
	dim1_ = m;
	dim2_ = n;
	Vector_hpc<T>::p_ = new T[m*n];
	Vector_hpc<T>::dim_ = m*n;
	zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_hpc<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc(matrix_manner_t type, integer_t m, integer_t n)
{
    if (m<=0 || n<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = dim2_ = 0;
        std::cerr << "Error: bad value in Matrix_hpc constructor " << std::endl;
        return;
    }
    if (type == ROW_MAJOR) {
        this->type = 0;
    } else {
        this->type = 1;
    }
    dim1_ = m;
    dim2_ = n;
    Vector_hpc<T>::p_ = new T[m*n];
    Vector_hpc<T>::dim_ = m*n;
    zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_hpc<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::¡«Matrix_hpc()
{
    std::cout << "Matrix_hpc destructor"<< std::endl;
}

template <typename T>
Matrix_hpc<T>& Matrix_hpc<T>::newsize(integer_t m, integer_t n)
{
    if (m<=0 || n<=0) {
        return;
    }
    if (Vector_hpc<T>::p_) delete [] Vector_hpc<T>::p_;
    dim1_ = m;
    dim2_ = n;
    Vector_hpc<T>::p_ = new T[m*n];
    if (Vector_hpc<T>::p_ == NULL)
    {
    	std::cerr << "Error : NULL pointer in operator= Matrix_hpc newsize" << std::endl;
        exit(-1);
    }
    Vector_hpc<T>::dim_ = m*n;
    zero();
    nonzeroes = 0;

    return *this;
}

template <typename T>
Matrix_hpc<T>& Matrix_hpc<T>::operator=(const Matrix_hpc<T>& M)
{
    integer_t N = M.size();
    integer_t i;

    if (&M == this) {
        return *this;
    }

    // no need to test for overlap, since this region is new
    for (i =0; i< N; i++)       // careful not to use bcopy()
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];        // here, but double::operator= double.
    nonzeroes = M.nonzeroes;

    return *this;
}

template <typename T>
Matrix_hpc<T>& Matrix_hpc<T>::operator=(const T& m)
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
void Matrix_hpc<T>::add(const Matrix_hpc<T> &c1)
{
    integer_t dim = size();
    
    if (dim <= 0 || dim != c1.size()) {
        return;
    }

    for (integer_t i=0; i<dim; i++) {
        Vector_hpc<T>::p_[i] += c1.Vector_hpc<T>::p_[i];
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const Matrix_hpc<T> &M)
{
    integer_t N = M.size();
    integer_t i_max = M.dim1();
    integer_t j_max = M.dim2();

    s << "\nNumber of Rows = " << i_max << std::endl;
    s << "\nNumber of Cols = " << j_max << std::endl;
    s << std::endl;
    s.width(10);
    s << "    Row Index ";
    s.width(10);
    s << "    Col Index ";
    s.width(20);
    s << "    Value";
    s << std::endl;

    for (integer_t i=0; i< i_max; i++) {
        for (integer_t j=0; j< j_max; j++) {
            s.width(10);
            s << i <<"    ";
            s.width(10);
            s << j <<"    ";
            s.width(20);
            s << M(i,j);
            s << std::endl;;
        }
    }
    
    s << std::endl;

    return s;
}


