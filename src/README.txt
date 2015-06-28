ADT(abstract data structure):

This directory contains various ADT, which could be used in reactor calculations.
There are mainly two basic ADT.
Cmv.h               Vector_hpc class, vector for serial and distributed environment
Matrix_hpc.h        Matrix_hpc class, matrix for serial and distributed environment
DenseMatrix.h       DenseMatrix class, matrix for serial and distributed environment
SparseMatrix.h      Sparse matrix class, with CSR store manner and so on

Then some useful ADTs for nuclear reactor calculation are chose.
Dimscal.h
CrossSection.h
Flux.h
Stencil2D.h
Stencil3D.h

The developer could utilize these classes to construct new class. 
