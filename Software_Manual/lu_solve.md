# Solve a Square Linear System Using LU Decomposition and Forward Substitution

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           lu_solve

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c lu_solve.f90```

and can be added to a program using

```$ gfortran program.f90 lu_solve.o ``` 

**Description/Purpose:** This routine takes a decomposed a square coefficient matrix, ***A = L\*U*** (like that returned by the [lu_factor](./lu_factor.md) subroutine) and implements a forward substitution such that

<a href="https://www.codecogs.com/eqnedit.php?latex=L\mathbf{\vv{y}}&space;=&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L\mathbf{\vv{y}}&space;=&space;\mathbf{\vv{b}}" title="L\mathbf{\vv{y}} = \mathbf{\vv{b}}" /></a>

to find ***y***, followed by a backward substitution

<a href="https://www.codecogs.com/eqnedit.php?latex=U\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{y}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{y}}" title="U\mathbf{\vv{x}} = \mathbf{\vv{y}}" /></a>

to find ***x***, the solution to the system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=A\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{b}}" title="A\mathbf{\vv{x}} = \mathbf{\vv{b}}" /></a>.

**Input:** 

*n* : INTEGER - size of the square matrix *LU* and the right-hand side vector *b*

*LU* : REAL - the square, LU-decomposed coefficient matrix of size *n*

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
CALL lu_factor(A, n)
! `A` now stores its LU factorization.
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
SUBROUTINE lu_solve(LU, n, b, x)
	IMPLICIT NONE
	
	! Takes as inputs a square, LU-decomposed, coefficient matrix, `LU`
	! of size `n` x `n` and the right-hand side vector `b` to solve the
	! system of equations Ax=b for `x` using LU decomposition with
	! forward and backward substitution.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: b(1:n)
	REAL*8, INTENT(INOUT) :: LU(1:n, 1:n), x(1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: fsum, backsum
	REAL*8 :: y(1:n)
	
	! Implement the forward substitution algorithm on Ly = b to find
	! `y`.
	
	y(1) = b(1)/LU(1, 1)
	
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + LU(k, j)*y(j)
		END DO
		y(k) = (b(k) - fsum)
	END DO
	
	! Implement the backward substitution algorithm on Ux = y to find
	! `x`.
	
	x(n) = y(n)/LU(n, n)
	
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + LU(k, j)*x(j)
		END DO
		x(k) = (y(k) - backsum)/LU(k, k)
	END DO
	
END SUBROUTINE
```



