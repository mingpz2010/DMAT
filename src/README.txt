ADT(abstract data structure):

This directory contains various ADT, which could be used in reactor calculations.
There are mainly two basic ADT.
(1)
Cmv.h               Vector_hpc class, vector for serial and distributed environment
Matrix_hpc.h        Matrix_hpc class, matrix for serial and distributed environment
DenseMatrix.h       DenseMatrix class, matrix for serial and distributed environment
SparseMatrix.h      Sparse matrix class, with CSR store manner and so on

Then some useful ADTs for nuclear reactor calculation are chose.
(2)
Dimscal.h					3D abstract data type
Flux.h						4D abstract data type
CrossSection.h				5D abstract data type
LegendreScal.h				6D abstract data type
Stencil2D.h					2D stencil computations
Stencil3D.h					3D stencil computations

The developer could utilize these classes to construct new class. 
