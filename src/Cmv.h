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


class Vector_double
{
    protected:
        double *p_;
        integer_t dim_;
        int ref_;           // 0: own memory space; 1: point to another memory space
    public:
        /*::::::::::::::::::::::::::*/
        /* Constructors/Destructors */
        /*::::::::::::::::::::::::::*/
        Vector_double();
        Vector_double(integer_t);
        Vector_double(integer_t, const double&);
        Vector_double(double*, integer_t);
        Vector_double(const double*, integer_t);
        Vector_double(const Vector_double &);
        ~Vector_double();
        /*::::::::::::::::::::::::::::::::*/
        /*  Indices and access operations */
        /*::::::::::::::::::::::::::::::::*/
        double& operator()(integer_t i) {
            return p_[i];
        }
        const double& operator()(integer_t i) const {
            return p_[i];
        }
        double& operator[](integer_t i) {
            return p_[i];
        }
        const double& operator[](integer_t i) const {
            return p_[i];
        }

        inline integer_t size() const { return dim_;}
        inline integer_t dim() const { return dim_;}
        inline integer_t ref() const { return ref_; }
        inline int null() const {return dim_== 0;}
        //
        // Create a new *uninitalized* vector of size N
        Vector_double & newsize(integer_t);
        /*::::::::::::::*/
        /*  Assignment  */
        /*::::::::::::::*/
        Vector_double & operator=(const Vector_double&);
        Vector_double & operator=(const double&);
        friend Vector_double operator+(const Vector_double &c1, const Vector_double &c2);
        friend Vector_double operator+(const Vector_double &c1, double num);
        friend Vector_double operator+(double num, const Vector_double &c1);

        // common functions
        void add(const Vector_double &c1);
        void add(double *);
        void sub(const Vector_double &c1);
        void sub(double *);
        void mul(double num);
        void div(double num);
        double max();
        double min();
        double mean();

        // something related to Fortran
        void copyFortran(int ref, REAL8 *, INTEGER dim);

        friend std::ostream& operator<<(std::ostream &s, const Vector_double &A);
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

