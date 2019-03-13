# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           lu_solve

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c lu_solve.f90```

and can be added to a program using

```$ gfortran program.f90 lu_solve.o ``` 

**Description/Purpose:** This routine decomposes a square coefficient matrix, ***A***, into its LU components via the LU decomposition method, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;a_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;a_{21}&space;&a_{22}&space;&a_{23}&space;\\&space;a_{31}&space;&a_{32}&space;&a_{33}&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;\ell_{21}&space;&1&space;&0&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&1&space;\end{bmatrix}\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;0&space;&u_{22}&space;&u_{23}&space;\\&space;0&space;&0&space;&u_{33}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;a_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;a_{21}&space;&a_{22}&space;&a_{23}&space;\\&space;a_{31}&space;&a_{32}&space;&a_{33}&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;\ell_{21}&space;&1&space;&0&space;\\&space;\ell_{31}&space;&\ell_{32}&space;&1&space;\end{bmatrix}\begin{bmatrix}&space;u_{11}&space;&u_{12}&space;&u_{13}&space;\\&space;0&space;&u_{22}&space;&u_{23}&space;\\&space;0&space;&0&space;&u_{33}&space;\end{bmatrix}" title="\begin{bmatrix} a_{11} &a_{12} &a_{13} \\ a_{21} &a_{22} &a_{23} \\ a_{31} &a_{32} &a_{33} \end{bmatrix} = \begin{bmatrix} 1 &0 &0 \\ \ell_{21} &1 &0 \\ \ell_{31} &\ell_{32} &1 \end{bmatrix}\begin{bmatrix} u_{11} &u_{12} &u_{13} \\ 0 &u_{22} &u_{23} \\ 0 &0 &u_{33} \end{bmatrix}" /></a>

and implements a forward substitution such that

<a href="https://www.codecogs.com/eqnedit.php?latex=L\mathbf{\vv{y}}&space;=&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L\mathbf{\vv{y}}&space;=&space;\mathbf{\vv{b}}" title="L\mathbf{\vv{y}} = \mathbf{\vv{b}}" /></a>

to find ***y***, followed by a backward substitution

<a href="https://www.codecogs.com/eqnedit.php?latex=U\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{y}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{y}}" title="U\mathbf{\vv{x}} = \mathbf{\vv{y}}" /></a>

to find ***x***, the solution to the system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=A\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{b}}" title="A\mathbf{\vv{x}} = \mathbf{\vv{b}}" /></a>.

**Input:** 

*n* : INTEGER - size of the square matrix *A* and the right-hand side vector *b*

*A* : REAL - the square coefficient matrix of size *n*

*b* : REAL - the right-hand side vector

**Output:** 

*x* : REAL - the solution to the system of equations contained in *A* and *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:)

n = 3
ALLOCATE(A(1:n, 1:n), b(1:n), x(1:n))
A = RESHAPE((/1.D0, -1.D0, 3.D0, &
	& 1.D0, 1.D0, 0.D0, &
	& 3.D0, -2.D0, 1.D0/), (/n, n/), ORDER=(/2, 1/))
b = (/2.D0, 4.D0, 1.D0/)
CALL lu_solve(A, n, b, x)
DO i = 1, n
	WRITE(*,*) x(i)
END DO
```

The outputs from the above code:

```fortran
1.6153846153846154     
2.3846153846153846     
0.92307692307692313 
```

**Implementation/Code:** The code for lu_solve can be seen below.

```fortran
SUBROUTINE lu_solve(A, n, b, x)
	IMPLICIT NONE
	
	! Takes as inputs a square coefficient matrix, `A` of size
	! `n` x `n` and the right-hand side vector `b` to solve the system
	! of equations Ax=b for `x` using LU decomposition with forward and
	! backward substitution.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: b(1:n)
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n), x(1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor, fsum, backsum
	REAL*8 :: y(1:n)
	
	! Find the LU decomposition of the coefficient matrix `A`.
	DO k = 1, n - 1
		
		DO i = k + 1, n
			
			factor = A(i, k)/A(k, k)
			
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
			
			A(i, k) = factor
		END DO
	END DO
	
	! Implement the forward substitution algorithm on Ly = b to find
	! `y`.
	
	y(1) = b(1)/A(1, 1)
	
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + A(k, j)*y(j)
		END DO
		y(k) = (b(k) - fsum)
	END DO
	
	! Implement the backward substitution algorithm on Ux = y to find
	! `x`.
	
	x(n) = y(n)/A(n, n)
	
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + A(k, j)*x(j)
		END DO
		x(k) = (y(k) - backsum)/A(k, k)
	END DO
	
END SUBROUTINE
```



