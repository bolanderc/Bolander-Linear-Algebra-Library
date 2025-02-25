# Solve a Square Linear (Lower-Triangular) System Using Forward Substitution

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           forwardsub

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c forwardsub.f90```

and can be added to a program using

```$ gfortran program.f90 forwardsub.o ``` 

**Description/Purpose:** This routine computes the solution of a square, linear system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

where the coefficient matrix, ***A*** is a lower-triangular matrix, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;a_{11}&space;&0&space;&0&space;\\&space;a_{21}&space;&a_{22}&space;&0&space;\\&space;a_{31}&space;&a_{32}&space;&a_{33}&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;b_1\\&space;b_2\\&space;b_3&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;a_{11}&space;&0&space;&0&space;\\&space;a_{21}&space;&a_{22}&space;&0&space;\\&space;a_{31}&space;&a_{32}&space;&a_{33}&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;b_1\\&space;b_2\\&space;b_3&space;\end{bmatrix}" title="\begin{bmatrix} a_{11} &0 &0 \\ a_{21} &a_{22} &0 \\ a_{31} &a_{32} &a_{33} \end{bmatrix}\begin{bmatrix} x_1\\ x_2\\ x_3 \end{bmatrix} = \begin{bmatrix} b_1\\ b_2\\ b_3 \end{bmatrix}" /></a>

using the forward substitution method, which has the algorithm

for *k* = 1 : *n*

<a href="https://www.codecogs.com/eqnedit.php?latex=x_k&space;=&space;\frac{b_k-\sum_{j=1}^{k-1}a_{kj}x_j}{a_{kk}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_k&space;=&space;\frac{b_k-\sum_{j=1}^{k-1}a_{kj}x_j}{a_{kk}}" title="x_k = \frac{b_k-\sum_{j=1}^{k-1}a_{kj}x_j}{a_{kk}}" /></a>

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix *A*, and the length of *b* and *x*

*A* : REAL - square, lower-triangular matrix of size *n* x *n*

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
A = RESHAPE((/3.D0, 0.D0, 0.D0, 0.D0, &
			& -1.D0, 6.D0, 0.D0, 0.D0, &
			& 3.D0, 2.D0, -16.D0, 0.D0, &
			& 1.D0, 1.D0, 1.D0, 1.D0/), (/n, n/), ORDER=(/2, 1/))
b = (/4.D0, 10.D0, 32.D0, 20.D0/)
CALL forwardsub(n, A, b, x)
DO i = 1, n
	WRITE(*,*) x(i)
END DO
```

The outputs from the above code:

```fortran
   1.3333333333333333     
   1.8888888888888891     
  -1.5138888888888888     
   18.291666666666668  
```

**Implementation/Code:** The code for forwardsub can be seen below.

```fortran
SUBROUTINE forwardsub(n, A, b, x)
	IMPLICIT NONE
	
	! Takes as input the size of the square matrix, `n`, the lower
	! triangular matrix, `A`, and the right-hand-side vector `b`.
	! Outputs the solution vector `x`.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Initialize decrement variables and a variable to compute the sum
	! of previous solutions integrated into the algorithm.
	INTEGER k, j
	REAL*8 fsum
	
	! Calculate the first value in the solution vector `x`.
	x(1) = b(1)/A(1, 1)
	
	! Loop through the remaining rows in `x` to calculate the solution
	! using the forward substitution algorithm.
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + A(k, j)*x(j)
		END DO
		x(k) = (b(k) - fsum)/A(k, k)
	END DO
END SUBROUTINE
```



