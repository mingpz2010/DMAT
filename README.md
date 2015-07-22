# DMAT

@(Matrix)[Mutiple Dimensional Array|Cache|Parallel Computing]

  It is a high efficiency data structure on vector, matrix and nuclear reactor physical quantities. It will be used to store physical variables and do operations in reactor calculation.

-------------------

[TOC]

##Abstract data type

> The reason for introducing the notion of abstract data types was to allow interchangeable software modules. You cannot have interchangeable modules unless these modules share similar complexity behavior. If I replace one module with another module with the same functional behavior but with different complexity tradeoffs, the user of this code will be unpleasantly surprised. I could tell him anything I like about data abstraction, and he still would not want to use the code. Complexity assertions have to be part of the interface.    ¡ª¡ª [Wikipedia](https://en.wikipedia.org/wiki/Abstract_data_type)

  As the vector and matrix operation are the basic of several scientific computing softwares. The C++ template of vector, matrix and multiple dimensional array will be researched. `MV++` s a small, efficient, set of concrete vector and simple matrix classes for numerical computing. It is not intended as a general vector container class, but rather designed specifically for optimized numerical computations on RISC and pipelined architectures. It is one step above a C/C++ array. `BLAS` are routines that provide standard building blocks for performing basic vector and matrix operations. The Level 1 BLAS perform scalar, vector and vector-vector operations, the Level 2 BLAS perform matrix-vector operations, and the Level 3 BLAS perform matrix-matrix operations. Because the BLAS are efficient, portable, and widely available, they are commonly used in the development of high quality linear algebra software, LAPACK for example. 

### Code Example
``` C++
template<typename T> class Dimscal;
template<typename T> std::ostream& operator<<(std::ostream &s, const Dimscal<T> &M);

template <typename T>
class Dimscal : public Vector_hpc<T>
{
protected:
    integer_t dim1_;
    integer_t dim2_;
    integer_t dim3_;
    integer_t w1, w2;
public:
    Dimscal() : Vector_hpc<T>() { dim1_= dim2_ = dim3_ = 0; w1 = w2 = 0; }
    Dimscal(integer_t, integer_t, integer_t);
    inline const T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) const {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }
    inline T& operator() (integer_t dim1, integer_t dim2,
            integer_t dim3) {
        return Vector_hpc<T>::p_[dim1*w1+dim2*w2+dim3];
    }
    inline T *ptr() { return Vector_hpc<T>::p_; }
    inline integer_t dim1() { return dim1_; }
    inline integer_t dim2() { return dim2_; }
    inline integer_t dim3() { return dim3_; }
    inline integer_t dim() { return Vector_hpc<T>::dim_; }
    inline void mul(T num) {
        for (integer_t i=0; i<Vector_hpc<T>::dim_; i++) {
            Vector_hpc<T>::p_[i] *= num;
        }
    }

    Dimscal<T> & operator=(const Dimscal<T>&);

    // common functions
    T maxPos(integer_t& dim1, integer_t& dim2, integer_t& dim3);
    void blas_op(T a, T b, T c);

    friend std::ostream& operator<< <>(std::ostream &s, const Dimscal<T> &M);
};

#include "Dimscal.tpp"
```
### Computing Formulation

Multiple Dimensional Array $B+c_{1}\times i_{1}+c_{2}\times i_{2}+...+c_{k}\times i_{k}$. And sometimes physical state variables in nuclear calculation are five or six dimensional array, such as neutron energy current values(Legendre polynomial expansion)£º

$$	\sum_{k=0}^\infty P_{n}(x)t^{n}=\frac{1}{\sqrt{1-2xt+t^{2}}} $$


### System Arichitecture
```flow
st=>start: Start
e=>end
op=>operation: My Operation
cond=>condition: Yes or No?

st->op->cond
cond(yes)->e
cond(no)->op
```



### Development Calendar
2015.5.1   Analysis MV++ NIST national library.
2015.6.2   Accomplish BLAS routines analysis and research some cache optimizations.
2015.7.1   Implement the basic library and use it in one 3D reactor core software.
2015.7.22  The BLAS routine was cutted and implemented in DMAT library.



## Feedback

- MAIL£º<mingpz@mail.ustc.edu.cn>

---------
Thank you for your attention to this library. Please download and utilize it in a appropriate style as your requirement.


