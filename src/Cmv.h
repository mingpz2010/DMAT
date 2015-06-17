/*
 * Cmv.h
 *
 *  Created on: 2015.6.11
 *      Author: PingzhouMing
 */

#ifndef CMV_H_
#define CMV_H_

#ifdef integer_t
#undef integer_t
#endif

typedef long integer_t;

//
// Fortran changed to C/C++ basic data type
//
//      C++ type        name          Fortran type
// ----------------------------------------------
typedef unsigned char   BYTE;       // byte
typedef short int       INTEGER2;   // integer*2
typedef long            INTEGER;    // integer
typedef int             LOGICAL;    // logical
typedef float           REAL;       // real
typedef double          REAL8;      // real*8
typedef long double     REAL16;     // real*16

// support float, double, int and long basic data type (generic programming)
template <typename T>
class Vector_hpc
{
    protected:
        T *p_;
        integer_t dim_;
        int ref_;           // 0: own memory space; 1: point to another memory space
    public:
        /*::::::::::::::::::::::::::*/
        /* Constructors/Destructors */
        /*::::::::::::::::::::::::::*/
        Vector_hpc() { p_ = NULL; dim_ = 0; }
        Vector_hpc(integer_t);
        Vector_hpc(integer_t, const T&);
        Vector_hpc(double*, integer_t);
        Vector_hpc(const double*, integer_t);
        Vector_hpc(const Vector_hpc &);
        ~Vector_hpc();
        /*::::::::::::::::::::::::::::::::*/
        /*  Indices and access operations */
        /*::::::::::::::::::::::::::::::::*/
        T& operator()(integer_t i) {
            return p_[i];
        }
        const T& operator()(integer_t i) const {
            return p_[i];
        }
        T& operator[](integer_t i) {
            return p_[i];
        }
        const T& operator[](integer_t i) const {
            return p_[i];
        }

        inline integer_t size() const { return dim_;}
        inline integer_t dim() const { return dim_;}
        inline integer_t ref() const { return ref_; }
        inline int null() const {return dim_== 0;}
        //
        // Create a new *uninitalized* vector of size N
        Vector_hpc & newsize(integer_t);
        /*::::::::::::::*/
        /*  Assignment  */
        /*::::::::::::::*/
        Vector_hpc & operator=(const Vector_hpc&);
        Vector_hpc & operator=(const T&);
        friend Vector_hpc operator+(const Vector_hpc &c1, const Vector_hpc &c2);
        friend Vector_hpc operator+(const Vector_hpc &c1, T num);
        friend Vector_hpc operator+(T num, const Vector_hpc &c1);

        // common functions
        void add(const Vector_hpc &c1);
        void add(T *);
        void sub(const Vector_hpc &c1);
        void sub(T *);
        void mul(T num);
        void div(T num);
        double max();
        double min();
        double mean();

        // something related to Fortran
        void copyFortran(int ref, T *, INTEGER dim);

        friend std::ostream& operator<<(std::ostream &s, const Vector_hpc &A);
};

class Face_current : public Vector_double
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t dim4_;
    integer_t w1, w2, w3;
public:
    Face_current() : Vector_double() { dim1_ = dim2_ = dim3_ = dim4_ = 0; }
    Face_current(integer_t dim1, integer_t dim2, integer_t dim3, integer_t dim4);
    inline const double& Face_current::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3, integer_t dim4) const {
        return p_[dim1*w1+dim2*w2+dim3*w3+dim4];
    }
    inline double& Face_current::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3, integer_t dim4) {
        return p_[dim1*w1+dim2*w2+dim3*w3+dim4];
    }
    ~Face_current();
};

class Dimensional_scal : public Vector_double
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t w1, w2;
public:
    Dimensional_scal() : Vector_double() { dim1_= dim2_ = dim3_ = 0; }
    Dimensional_scal(integer_t dim1, integer_t dim2, integer_t dim3);
    inline const double& Dimensional_scal::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) const {
        return p_[dim1*w1+dim2*w2+dim3];
    }
    inline double& Dimensional_scal::operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) {
        return p_[dim1*w1+dim2*w2+dim3];
    }
    ~Dimensional_scal();

    // common functions
    double maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3);

    // something related to Fortran
    void copyFortran(int ref, REAL8 *, INTEGER dim1, INTEGER dim2, INTEGER dim3);
};


#endif /* CMV_H_ */

