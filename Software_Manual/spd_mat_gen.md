# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           spd_mat_gen

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c spd_mat_gen.f90```

and can be added to a program using

```$ gfortran program.f90 spd_mat_gen.o ``` 

**Description/Purpose:** This routine takes a square matrix and turns it into a symmetric, positive definite (SPD) matrix. A symmetric, positive definite matrix ***A*** can be found by filling a matrix ***C***, of the same size as ***A***, with random numbers. As long as ***C*** is not singular, ***A*** can be found by

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{C}}^T\mathbf{\vv{C}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{C}}^T\mathbf{\vv{C}}" title="\mathbf{\vv{A}} = \mathbf{\vv{C}}^T\mathbf{\vv{C}}" /></a>

A test to ensure that ***C*** is not singular is to run it through the Cholesky factorization method, which is outlined in the [cholesky_factor](./cholesky_factor.md) subroutine. This is performed in the subroutine.

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix

**Output:** 

*A* : REAL - a symmetric, positive definite array of size (n, n)

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i, LDA, LWMAX, INFO, LWORK
REAL*8, ALLOCATABLE :: A(:, :), W(:), WORK(:)

n = 3
LWMAX = 1000
LDA = 3
ALLOCATE(A(1:n, 1:n), W(1:n), WORK(1:LWMAX))
CALL spd_mat_gen(A, n)
DO i = 1, n
	WRITE(*,*) A(i, :)
END DO

! Find eigenvalues to see if this is really a SDM.
LWORK = -1
CALL DSYEV( 'N', 'U', N, A, LDA, W, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DSYEV( 'N', 'U', N, A, LDA, W, WORK, LWORK, INFO )

IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
WRITE(*,*) W
```

The outputs from the above code for several iterations:

```fortran
 A (SPD)
 --------------
  0.95913044795548175        1.0770792618925378       0.74378929487883871     
   1.0770792618925378        1.5533761602066278       0.90661535805392335     
  0.74378929487883871       0.90661535805392335        1.1326963806643067     
 Eigenvalues (to prove this is an SPD matrix)
 --------------
 0.13425525842489663       0.43164312395922544        3.0793046064422933  
  
 A (SPD)
 --------------
   1.1311766249561668       0.74513887812700086       0.82047235318646705     
  0.74513887812700086       0.58482114584560896       0.44055505581139276     
  0.82047235318646705       0.44055505581139276       0.74440500606864424     
 Eigenvalues (to prove this is an SPD matrix)
 --------------
9.2733773949918285E-003  0.22142261442256186        2.2297067850528656

 A (SPD)
 --------------
1.1232444860161654   1.0002849940748892  0.79601097869314197   1.2138409975251179     
1.0002849940748892   1.1800569306915463  0.52850657698680170   1.2585964158335812     
0.79601097869314197  0.52850657698680170  0.91304479758715851  0.79921994769253413     
1.2138409975251179   1.2585964158335812  0.79921994769253413   1.4695759719362007     
 Eigenvalues (to prove this is an SPD matrix)
 --------------
2.2277439051183831E-002 6.1118892611140660E-002 0.53050555974032199 4.0720202948284241 
```

An additional, large example can be seen [here](../HW5Task4Report.md).

**Implementation/Code:** The code for spd_mat_gen can be seen below.

```fortran
SUBROUTINE spd_mat_gen(A, n)
	IMPLICIT NONE
	
	! Takes as input the size `n` of the square matrix `A` and outputs
	! a symmetric, positive definite matrix in `A`.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: A(1:n, 1:n)
	
	REAL*8 :: C(1:n, 1:n)
	INTEGER :: error
	
	error = 1
	
	! Make sure it is a proper spd matrix.
	DO WHILE (error == 1)
		! Creates a matrix `C` filled with random numbers as well as its
		! transpose.
		CALL RANDOM_NUMBER(C)
		
		! Since a symmetric, positive definite matrix can be defined as a
		! matrix's transpose multiplied by original matrix, (i.e. A = CT*C)
		! `A` is determined as a product of these two matrices. Note that
		! this will only fail to produce an SPD matrix if C is singular.
		CALL mat_prod(TRANSPOSE(C), C, n, n, n, A)
		
		! Check to see if the matrix is really symmetric positive definite
		! by running it through the `cholesky_factor` subroutine.
		CALL cholesky_factor(A, n, error)
	END DO

END SUBROUTINE
```
