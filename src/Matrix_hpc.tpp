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
    property = M_NORMAL;
    Vector_hpc<T>::p_ = NULL;
    Vector_hpc<T>::dim_ = 0;
	dim1_ = 0;
	nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc(integer_t n)
{
	if (m<=0 || n<=0) {
	    property = M_NORMAL;
	    Vector_hpc<T>::p_ = NULL;
	    Vector_hpc<T>::dim_ = 0;
	    dim1_ = 0;
	    std::cerr << "Error: bad value in Matrix_hpc constructor " << std::endl;
	    return;
	}
	type = 0;
	property = M_NORMAL;
	dim1_ = n;
	Vector_hpc<T>::p_ = new T[n*n];
	Vector_hpc<T>::dim_ = n*n;
	zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_hpc<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::Matrix_hpc(matrix_manner_t type, integer_t n)
{
    if (m<=0 || n<=0) {
        Vector_hpc<T>::p_ = NULL;
        Vector_hpc<T>::dim_ = 0;
        dim1_ = 0;
        std::cerr << "Error: bad value in Matrix_hpc constructor " << std::endl;
        return;
    }
    if (type == ROW_MAJOR) {
        this->type = 0;
    } else {
        this->type = 1;
    }
    property = M_NORMAL;
    dim1_ = n;
    Vector_hpc<T>::p_ = new T[n*n];
    Vector_hpc<T>::dim_ = n*n;
    zero();
    if (Vector_hpc<T>::p_ == NULL) {
        std::cerr << "Error: NULL pointer in Matrix_hpc<T> constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    nonzeroes = 0;
}

template <typename T>
Matrix_hpc<T>::~Matrix_hpc()
{

}

template <typename T>
Matrix_hpc<T>& Matrix_hpc<T>::newsize(integer_t n)
{
    if (n<=0) {
        return *this;
    }
    if (Vector_hpc<T>::p_) delete [] Vector_hpc<T>::p_;
    dim1_ = n;
    Vector_hpc<T>::p_ = new T[n*n];
    if (Vector_hpc<T>::p_ == NULL)
    {
    	std::cerr << "Error : NULL pointer in operator= Matrix_hpc newsize" << std::endl;
        exit(-1);
    }
    Vector_hpc<T>::dim_ = n*n;
    zero();
    nonzeroes = 0;

    return *this;
}

//
//  Here is the BLAS-cutted interfaces
//
template <typename T>
void Matrix_hpc<T>::blas_op(T a, T b, T c)
{
    if (c == 0)  return;

    for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
        Vector_hpc<T>::p_[i] *= a;
        Vector_hpc<T>::p_[i] += b;
        Vector_hpc<T>::p_[i] /= c;
    }
}

template <typename T>
void Matrix_hpc<T>::dswap(const Matrix_hpc<T>& M)
{
    if (dim_ <= 0 || dim_ != M.dim())  return;

    T tmp;
    for ( integer_t i=0; i<dim_; i++) {
        T = M.Vector_hpc<T>::p_[i];
        M.Vector_hpc<T>::p_[i] = Vector_hpc<T>::p_[i];
        Vector_hpc<T>::p_[i] = T;
    }
}

template <typename T>
void Matrix_hpc<T>::dscal(T a)
{
    // unroll loops to depth of length 8
    integer_t N = size();
    integer_t Nminus8 = N-8;
    integer_t i;

    for (i=0; i<Nminus8; )
    {
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
        Vector_hpc<T>::p_[i] *= a; ++i;
    }

    for (; i<N; Vector_hpc<T>::p_[i] *= a, ++i);   // finish off last piece...
}

template <typename T>
void Matrix_hpc<T>::dcopy(const Matrix_hpc<T>& M)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i];  ++i;
    }

    for (; i<N; Vector_hpc<T>::p_[i] = M.Vector_hpc<T>::p_[i], ++i);  // finish off last piece...
}

template <typename T>
void Matrix_hpc<T>::daxpy(T a, const Matrix_hpc<T>& M)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        Vector_hpc<T>::p_[i] += a*M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] += a*M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] += a*M.Vector_hpc<T>::p_[i];  ++i;
        Vector_hpc<T>::p_[i] += a*M.Vector_hpc<T>::p_[i];  ++i;
    }

    for (; i<N; ++i) {
        // finish off last piece...
        Vector_hpc<T>::p_[i] += a*M.Vector_hpc<T>::p_[i];
    }
}

template <typename T>
T Matrix_hpc<T>::dnrm2()
{
    T ans = 0.;
    for ( integer_t i=0; i<dim_; i++) {
        ans += (Vector_hpc<T>::p_[i]*Vector_hpc<T>::p_[i]);
    }

    return sqrt(ans);
}

template <typename T>
T Matrix_hpc<T>::dasum()
{
    T ans = 0.;
    for ( integer_t i=0; i<dim_; i++) {
        ans += fabs(Vector_hpc<T>::p_[i]);
    }

    return ans;
}

template <typename T>
void Matrix_hpc<T>::dgemv(Vector_hpc<T>& v)
{
    if (dim1_ != v.dim()) {
        std::cerr << "DGEMV dimension is not match!" << std::endl;
        return;
    }

    Vector_hpc<T> v_tmp(v);
    T ans;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<dim1_; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v_tmp[j];
            }
            v[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<dim1_; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v_tmp[j];
            }
            v[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dgemv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2)
{
    if (dim1_ != v1.dim() || dim1_ != v2.dim()) {
        std::cerr << "DGEMV dimension is not match!" << std::endl;
        return;
    }
    T ans;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<dim1_; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v1[j];
            }
            v2[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<dim1_; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v1[j];
            }
            v2[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dgbmv(Vector_hpc<T>& v)
{
    if (property != M_BANDED && property != M_SYMMETRIC_BANDED && property != M_TRIANGULAR) {
        return;
    }
    if (dim1_ != v.dim()) {
        std::cerr << "DGBMV dimension is not match!" << std::endl;
        return;
    }

    Vector_hpc<T> v_tmp(v);
    T ans;
    integer_t index;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            index = ((i-ksub) < 0)? 0: i-ksub;
            for (integer_t j=index; j<i; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v_tmp[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v_tmp[i];
            index = ((i+ksupper)>dim1_) ? dim1_:i+ksupper;
            for (integer_t j=i+1; j<index; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v_tmp[j];
            }
            v[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            index = ((i-ksub) < 0)? 0: i-ksub;
            for (integer_t j=index; j<i; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v_tmp[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v_tmp[i];
            index = ((i+ksupper)>dim1_) ? dim1_:i+ksupper;
            for (integer_t j=i+1; j<index; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v_tmp[j];
            }
            v[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dgbmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2)
{
    if (property != M_BANDED && property != M_SYMMETRIC_BANDED && property != M_TRIANGULAR) {
        return;
    }

    if (dim1_ != v1.dim() || dim1_ != v2.dim()) {
        std::cerr << "DGBMV dimension is not match!" << std::endl;
        return;
    }
    T ans;
    integer_t index;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            index = ((i-ksub) < 0)? 0: i-ksub;
            for (integer_t j=index; j<i; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v1[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v1[i];
            index = ((i+ksupper)>dim1_) ? dim1_:i+ksupper;
            for (integer_t j=i+1; j<index; j++) {
                ans += Vector_hpc<T>::p_[i*dim1_ + j]*v1[j];
            }
            v2[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            index = ((i-ksub) < 0)? 0: i-ksub;
            for (integer_t j=index; j<i; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v1[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v1[i];
            index = ((i+ksupper)>dim1_) ? dim1_:i+ksupper;
            for (integer_t j=i+1; j<index; j++) {
                ans += Vector_hpc<T>::p_[i + j*dim1_]*v1[j];
            }
            v2[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dsymv(Vector_hpc<T>& v)
{
    if (property != M_SYMMETRIC) {
        return;
    }
    if (dim1_ != v.dim()) {
        std::cerr << "DSYMV dimension is not match!" << std::endl;
        return;
    }

    Vector_hpc<T> v_tmp(v);
    T ans;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<i; j++) {
                ans += 2*Vector_hpc<T>::p_[i*dim1_ + j]*v_tmp[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v_tmp[i];
            v[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<i; j++) {
                ans += 2*Vector_hpc<T>::p_[i + j*dim1_]*v_tmp[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v_tmp[i];
            v[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dsymv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2)
{
    if (property != M_SYMMETRIC) {
        return;
    }

    if (dim1_ != v1.dim() || dim1_ != v2.dim()) {
        std::cerr << "DSYMV dimension is not match!" << std::endl;
        return;
    }
    T ans;
    if (type == 0) {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<i; j++) {
                ans += 2*Vector_hpc<T>::p_[i*dim1_ + j]*v1[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v1[i];
            v2[i] = ans;
        }
    } else {
        for (integer_t i=0; i<dim1_; i++) {
            ans = 0.;
            for (integer_t j=0; j<=i; j++) {
                ans += 2*Vector_hpc<T>::p_[i + j*dim1_]*v1[j];
            }
            ans += Vector_hpc<T>::p_[i*dim1_ + i]*v1[i];
            v2[i] = ans;
        }
    }
}

template <typename T>
void Matrix_hpc<T>::dtrmv(Vector_hpc<T>& v)
{
    if (property != M_TRIANGULAR) {
        return;
    }
    if (dim1_ != v.dim()) {
        std::cerr << "DTRMV dimension is not match!" << std::endl;
        return;
    }

    dgbmv(v);
}

template <typename T>
void Matrix_hpc<T>::dtrmv(const Vector_hpc<T>& v1, Vector_hpc<T>& v2)
{
    if (property != M_TRIANGULAR) {
        return;
    }

    if (dim1_ != v1.dim() || dim1_ != v2.dim()) {
        std::cerr << "DTRMV dimension is not match!" << std::endl;
        return;
    }

    dgbmv(v1, v2);
}

template <typename T>
void Matrix_hpc<T>::dtrsv(Vector_hpc<T>& x, const Vector_hpc<T>& b)
{
    if (property != M_TRIANGULAR) {
        return;
    }
    if (dim1_ != x.dim() || x.dim() != b.dim()) {
        std::cerr << "DTRSV dimension is not match!" << std::endl;
        return;
    }
    if (dim1_ == 1) {
        x[0] = b[0]/Vector_hpc<T>::p_[0];
        return;
    }
    if (dim1_ == 2) {
        x[1] = b[1]/Vector_hpc<T>::p_[3];
        x[0] = (b[0]-Vector_hpc<T>::p_[1]*x[1])/Vector_hpc<T>::p_[0];
        return;
    }

    {
        // GONO!
    }
}

template <typename T>
void Matrix_hpc<T>::dgemm(Matrix_hpc<T>& m)
{

}

template <typename T>
void Matrix_hpc<T>::dgemm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2)
{

}

template <typename T>
void Matrix_hpc<T>::dsymm(Matrix_hpc<T>& m)
{

}

template <typename T>
void Matrix_hpc<T>::dsymm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2)
{

}

template <typename T>
void Matrix_hpc<T>::dtrmm(Matrix_hpc<T>& m)
{

}

template <typename T>
void Matrix_hpc<T>::dtrmm(const Matrix_hpc<T>& m1, Matrix_hpc<T>& m2)
{

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
    integer_t j_max = M.dim1();

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


