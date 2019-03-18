# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           lu_factor

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c lu_factor.f90```

and can be added to a program using

```$ gfortran program.f90 lu_factor.o ``` 

**Description/Purpose:** This routine decomposes a square coefficient matrix, ***A*** into its LU factorization, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{L}}\mathbf{\vv{U}}&space;=&space;\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;\ell_{21}&space;&1&space;&0&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&1&space;\end{bmatrix}&space;\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;0&space;&u_{22}&space;&u_{23}&space;\\&space;0&space;&0&space;&u_{33}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{L}}\mathbf{\vv{U}}&space;=&space;\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;\ell_{21}&space;&1&space;&0&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&1&space;\end{bmatrix}&space;\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;0&space;&u_{22}&space;&u_{23}&space;\\&space;0&space;&0&space;&u_{33}&space;\end{bmatrix}" title="\mathbf{\vv{A}} = \mathbf{\vv{L}}\mathbf{\vv{U}} = \begin{bmatrix} 1 &0 &0 \\ \ell_{21} &1 &0 \\ \ell_{31} &\ell_{32} &1 \end{bmatrix} \begin{bmatrix} u_{11} &u_{12} &u_{13} \\ 0 &u_{22} &u_{23} \\ 0 &0 &u_{33} \end{bmatrix}" /></a>

and stores it as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{LU}}&space;=&space;\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;\ell_{21}&space;&u_{22}&space;&u_{23}&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&u_{33}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{LU}}&space;=&space;\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;\ell_{21}&space;&u_{22}&space;&u_{23}&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&u_{33}&space;\end{bmatrix}" title="\mathbf{\vv{LU}} = \begin{bmatrix} u_{11} &u_{12} &u_{13} \\ \ell_{21} &u_{22} &u_{23} \\ \ell_{31} &\ell_{32} &u_{33} \end{bmatrix}" /></a>

This can be used as part of the [lu_solve](./lu_solve.md) subroutine to solve a linear system of equations using the LU factorization method.

**Input:** 

*n* : INTEGER - size of the square matrix *LU* and the right-hand side vector *b*

*A* : REAL - the square coefficient matrix of size *n*

*b* : REAL - the right-hand side vector

**Output:** 

*x* : REAL - the solution to the system of equations contained in *A* and *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: A(:, :)

n = 3
ALLOCATE(A(1:n, 1:n))
A = RESHAPE((/1.D0, -1.D0, 3.D0, &
	& 1.D0, 1.D0, 0.D0, &
	& 3.D0, -2.D0, 1.D0/), (/n, n/), ORDER=(/2, 1/))
CALL lu_factor(A, n)
DO i = 1, n
	WRITE(*,*) A(i, :)
END DO
```

The outputs from the above code:

```fortran
   1.0000000000000000       -1.0000000000000000        3.0000000000000000     
   1.0000000000000000        2.0000000000000000       -3.0000000000000000     
   3.0000000000000000       0.50000000000000000       -6.5000000000000000 
```

**Implementation/Code:** The code for lu_factor can be seen below.

```fortran
SUBROUTINE lu_factor(A, n)
	IMPLICIT NONE
	
	! Takes as inputs a square coefficient matrix, `A` of size
	! `n` x `n` and outputs the LU factorization of that matrix.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! Loops through all columns except for the last.
	DO k = 1, n - 1
		
		! Loops through the rows in the matrix `A`.
		DO i = k + 1, n
			
			! Calculates the factor to be used in both the Gaussian
			! Elimination algorithm used to find the U matrix and 
			! stored as the L matrix.
			factor = A(i, k)/A(k, k)
			
			! Performs the Gaussian Elimination and finds the U matrix
			! components.
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
			
			! Stores the factor used in the corresponding L matrix
			! component.
			A(i, k) = factor
		END DO
	END DO	
END SUBROUTINE
```



