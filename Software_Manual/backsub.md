# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           backsub

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c backsub.f90```

and can be added to a program using

```$ gfortran program.f90 backsub.o ``` 

**Description/Purpose:** This routine computes the solution of a square, linear system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

where the coefficient matrix, ***A*** is an upper-triangular matrix, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;a_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;0&space;&a_{22}&space;&a_{23}&space;\\&space;0&space;&0&space;&a_{33}&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;b_1\\&space;b_2\\&space;b_3&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;a_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;0&space;&a_{22}&space;&a_{23}&space;\\&space;0&space;&0&space;&a_{33}&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;b_1\\&space;b_2\\&space;b_3&space;\end{bmatrix}" title="\begin{bmatrix} a_{11} &a_{12} &a_{13} \\ 0 &a_{22} &a_{23} \\ 0 &0 &a_{33} \end{bmatrix}\begin{bmatrix} x_1\\ x_2\\ x_3 \end{bmatrix} = \begin{bmatrix} b_1\\ b_2\\ b_3 \end{bmatrix}" /></a>

using the backward substitution method, which has the algorithm

for *k* = *n* : -1 : 1

<a href="https://www.codecogs.com/eqnedit.php?latex=x_k&space;=&space;\frac{b_k-\sum_{j=k&plus;1}^na_{kj}x_j}{a_{kk}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_k&space;=&space;\frac{b_k-\sum_{j=k&plus;1}^na_{kj}x_j}{a_{kk}}" title="x_k = \frac{b_k-\sum_{j=k+1}^na_{kj}x_j}{a_{kk}}" /></a>

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix *A*, and the length of *b* and *x*

*A* : REAL - square, upper-triangular matrix of size *n* x *n*

*b* : REAL - arbitrary vector of length *n*

**Output:** 

*x* : REAL - the solution vector of length *n*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:)

n = 4
ALLOCATE(A(1:n, 1:n), b(1:n), x(1:n))
A = RESHAPE((/1.D0, 1.D0, 1.D0, 1.D0, &
			& 0.D0, -2.D0, -1.D0, -1.D0, &
			& 0.D0, 0.D0, 1.D0, -1.D0, &
			& 0.D0, 0.D0, 0.D0, -2.D0/), (/n, n/), ORDER=(/2, 1/))
b = (/4.D0, 3.D0, 2.D0, -7.D0/)
CALL backsub(n, A, b, x)
DO i = 1, n
	WRITE(*,*) x(i)
END DO
```

The outputs from the above code:

```fortran
   1.0000000000000000     
  -6.0000000000000000     
   5.5000000000000000     
   3.5000000000000000   
```

**Implementation/Code:** The code for backsub can be seen below.

```fortran
SUBROUTINE backsub(n, A, b, x)
	IMPLICIT NONE
	
	! Takes as an input the size of the square matrix `n`, the upper
	! triangular matrix `A`, and the right hand side vector `b`. Outputs
	! the solution vector `x`.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Initialize decrement variables and a variable to compute the sum
	! of previous solutions integrated into the algorithm.
	INTEGER :: k, j
	REAL*8 :: backsum
	
	! Calculate the last value in the solution vector `x`.
	x(n) = b(n)/A(n, n)
	
	! Loop through the remaining rows in `x` to calculate the solution
	! using the backward substitution algorithm.
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + A(k, j)*x(j)
		END DO
		x(k) = (b(k) - backsum)/A(k, k)
	END DO
END SUBROUTINE
```



