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
Vector_hpc<T>::Vector_hpc(integer_t n) : p_(new T[n]), dim_(n)
{
    if (p_ == NULL) {
        std::cerr << "Error: NULL pointer in Vector_hpc constructor " << std::endl;
        std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    zero();
}

template <typename T>
Vector_hpc<T>::Vector_hpc(integer_t n, const T& v) : p_(new T[n]), dim_(n)
{
    if (p_ == NULL)
    {
    	std::cerr << "Error: NULL pointer in Vector_hpc constructor " << std::endl;
    	std::cerr << "       Most likely out of memory... " << std::endl;
        exit(-1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = v;
}

template <typename T>
Vector_hpc<T>::Vector_hpc(T* d, integer_t n) : p_(new T[n]), dim_(n)
{
    if (p_ == NULL)
    {
    	std::cerr << "Error: Null pointer in Vector_hpc(double*, int) " << std::endl;
        exit(1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = d[i];

}

template <typename T>
Vector_hpc<T>::Vector_hpc(const T* d, integer_t n) : p_(new T[n]), dim_(n)
{
    if (p_ == NULL)
    {
    	std::cerr << "Error: Null pointer in Vector_hpc(double*, int) " << std::endl;
        exit(-1);
    }
    for (integer_t i=0; i<n; i++)
        p_[i] = d[i];

}

template <typename T>
Vector_hpc<T>::Vector_hpc(const Vector_hpc<T> & m) : p_(new T[m.dim_]), dim_(m.dim_)
{
    if (p_ == NULL)
    {
    	std::cerr << "Error:  Null pointer in Vector_hpc(const Marray_1D_double&); " << std::endl;
        exit(-1);
    }

    integer_t N = m.dim_;

    for (integer_t i=0; i<N; i++)
        p_[i] = m.p_[i];
}

template <typename T>
Vector_hpc<T>::~Vector_hpc()
{
#ifdef DEBUG
    std::cout << "Vector_hpc destructor"<< std::endl;
#endif
    if (p_) delete [] p_;
}

template <typename T>
Vector_hpc<T>& Vector_hpc<T>::newsize(integer_t n)
{
    if (dim_ != n )                     // only delete and new if
    {                                   // the size of memory is really
        if (p_) delete [] p_;           // changing, otherwise just
        p_ = new T[n];              // copy in place.
        if (p_ == NULL)
        {
        	std::cerr << "Error : NULL pointer in operator= newsize" << std::endl;
            exit(-1);
        }
        dim_ = n;
        zero();
    }

    return *this;
}

template <typename T>
Vector_hpc<T>& Vector_hpc<T>::operator=(const T & m)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
    }

    for (; i<N; p_[i++] = m);   // finish off last piece...

    return *this;
}

template <typename T>
Vector_hpc<T>& Vector_hpc<T>::operator=(const Vector_hpc<T> & m)
{

    integer_t N = m.dim_;
    integer_t i;

    if (&m == this) {
    	return *this;
    }

    // no need to test for overlap, since this region is new
    for (i =0; i< N; i++)       // careful not to use bcopy()
        p_[i] = m.p_[i];        // here, but double::operator= double.

    return *this;
}

template <typename T>
Vector_hpc<T> operator+(const Vector_hpc<T> &c1, const Vector_hpc<T> &c2)
{
    Vector_hpc<T> c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] = c1[i] + c2[i];
    }

    return c;
}

template <typename T>
Vector_hpc<T> operator+(const Vector_hpc<T> &c1, T num)
{
    Vector_hpc<T> c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] += num;
    }

    return c;
}

template <typename T>
Vector_hpc<T> operator+(T num, const Vector_hpc<T> &c1)
{
    Vector_hpc<T> c(c1.dim_);

    for (integer_t i=0; i < c1.dim_; i++) {
        c[i] += num;
    }

    return c;
}

template <typename T>
void Vector_hpc<T>::add(const Vector_hpc<T> &c1)
{
    if (dim_ <= 0 || dim_ != c1.dim_) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] += c1[i];
    }
}

template <typename T>
void Vector_hpc<T>::add(T *m)
{
    if (dim_ <= 0) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] += m[i];
    }
}

template <typename T>
void Vector_hpc<T>::sub(const Vector_hpc &c1)
{
    if (dim_ <= 0 || dim_ != c1.dim_) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] -= c1[i];
    }
}

template <typename T>
void Vector_hpc<T>::sub(T *m)
{
    if (dim_ <= 0) {
        return;
    }

    for (integer_t i=0; i<dim_; i++) {
        p_[i] -= m[i];
    }
}

template <typename T>
void Vector_hpc<T>::mul(T num)
{
    if (dim_ <= 0 )  return;

    integer_t i;

    for (i=0; i<dim_; i++) {
        p_[i] *= num;
    }
}

template <typename T>
void Vector_hpc<T>::div(T num)
{
    if (dim_ <= 0 )  return;

    if (num == 0) {
    	std::cerr << "div 0 error, so we don't process it"<< std::endl;
        return;
    }

    integer_t i;

    for (i=0; i<dim_; i++) {
        p_[i] /= num;
    }
}


template <typename T>
void Vector_hpc<T>::fill(T num)
{
    // unroll loops to depth of length 8
    integer_t N = size();
    integer_t Nminus8 = N-8;
    integer_t i;

    for (i=0; i<Nminus8; )
    {
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
        p_[i++] = num;
    }

    for (; i<N; p_[i++] = num);   // finish off last piece...
}

template <typename T>
T Vector_hpc<T>::max()
{
    if (dim_ <= 0) { return 0.; }
    if (dim_ == 1) { return p_[0]; }

    T ans;
    integer_t i;
    ans = p_[0];
    for (i=1; i<dim_; i++) {
        if( ans < p_[i]) ans = p_[i];
    }
    return ans;
}

template <typename T>
T Vector_hpc<T>::min()
{
    if (dim_ <= 0) { return 0.; }
    if (dim_ == 1) { return p_[0]; }

    T ans;
    integer_t i;
    ans = p_[0];
    for (i=1; i<dim_; i++) {
        if( ans > p_[i]) ans = p_[i];
    }
    return ans;
}

template <typename T>
T Vector_hpc<T>::mean()
{
    if (dim_ <= 0)  return 0.;

    T ans;
    integer_t i;
    for ( ans = 0., i=0; i<dim_; i++) {
        ans += p_[i];
    }

    return ans/dim_;
}

//
//  Here is the BLAS-cutted interfaces
//
template <typename T>
void Vector_hpc<T>::dswap(const Vector_hpc<T>& M)
{
    if (dim_ <= 0 || dim_ != M.dim())  return;

    T tmp;
    for ( integer_t i=0; i<dim_; i++) {
        T = M[i];
        M[i] = p_[i];
        p_[i] = T;
    }
}

template <typename T>
void Vector_hpc<T>::dscal(T a)
{
    // unroll loops to depth of length 8
    integer_t N = size();
    integer_t Nminus8 = N-8;
    integer_t i;

    for (i=0; i<Nminus8; )
    {
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
        p_[i] *= a; ++i;
    }

    for (; i<N; p_[i] *= a, ++i);   // finish off last piece...
}

template <typename T>
void Vector_hpc<T>::dcopy(const Vector_hpc<T>& M)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        p_[i] = M[i];  ++i;
        p_[i] = M[i];  ++i;
        p_[i] = M[i];  ++i;
        p_[i] = M[i];  ++i;
    }

    for (; i<N; p_[i] = M[i], ++i);   // finish off last piece...
}

template <typename T>
void Vector_hpc<T>::daxpy(T a, const Vector_hpc<T>& M)
{
    // unroll loops to depth of length 4
    integer_t N = size();
    integer_t Nminus4 = N-4;
    integer_t i;

    for (i=0; i<Nminus4; )
    {
        p_[i] += a*M[i];  ++i;
        p_[i] += a*M[i];  ++i;
        p_[i] += a*M[i];  ++i;
        p_[i] += a*M[i];  ++i;
    }

    for (; i<N; p_[i] += a*M[i], ++i);   // finish off last piece...
}

template <typename T>
T Vector_hpc<T>::ddot(const Vector_hpc<T>& M)
{
    if (dim_ <= 0 || dim_ != M.dim())  return 0;

    T ans = 0.;
    for ( integer_t i=0; i<dim_; i++) {
        ans += p_[i] * M[i];
    }

    return ans;
}

template <typename T>
T Vector_hpc<T>::dnrm2()
{
    T ans = 0.;
    for ( integer_t i=0; i<dim_; i++) {
        ans += p_[i]*p_[i];
    }

    return sqrt(ans);
}

template <typename T>
T Vector_hpc<T>::dasum()
{
    T ans = 0.;
    for ( integer_t i=0; i<dim_; i++) {
        ans += fabs(p_[i]);
    }

    return ans;
}

template <typename T>
integer_t Vector_hpc<T>::idamax(T& ans)
{
    if (dim_ <= 0) {
        return -1;
    }

    ans = p_[0];
    integer_t index = 0;
    for (integer_t i=0; i<dim_; i++) {
        if (fabs(p_[i]) > ans) {
            ans = fabs(p_[i]);
            index = i;
        }
    }

    return index;
}

template <typename T>
void Vector_hpc<T>::copyFortran(int ref, T *from, INTEGER dim)
{
    ref_ = ref;

    if (p_ != NULL) {
        delete []p_;
    }
    dim_ = dim;
    if (ref_ == 0) {
        p_ = new T[dim_];
    } else {
        p_ = from;
    }
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const Vector_hpc<T>& V)
{
    integer_t N = V.size();

    for (integer_t i=0; i< N; i++)
        s << V(i) << " ";
    
    s << std::endl;

    return s;
}

template <typename T>
T max_NRM2(integer_t N, const Vector_hpc<T> &x1, const Vector_hpc<T> &x2)
{
    T ans;

    ans = fabs((x1(0) - x2(0))/(x1(0)));
    for (integer_t i=1; i<N; i++) {
        if (fabs((x1(i) - x2(i))/(x1(i))) > ans) {
            ans = fabs((x1(i) - x2(i))/(x1(i)));
        }
    }

    return ans;
}

