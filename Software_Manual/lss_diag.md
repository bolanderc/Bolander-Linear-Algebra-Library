# Solve a Square (Diagonal) Linear System

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           lss_diag

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c lss_diag.f90```

and can be added to a program using

```$ gfortran program.f90 lss_diag.o ``` 

**Description/Purpose:** This routine computes the solution of a square, linear system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

where the coefficient matrix, ***A*** is a diagonal matrix, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;0&space;&2&space;&0&space;\\&space;0&space;&0&space;&3&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;1\\&space;2\\&space;3&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;1&space;&0&space;&0&space;\\&space;0&space;&2&space;&0&space;\\&space;0&space;&0&space;&3&space;\end{bmatrix}\begin{bmatrix}&space;x_1\\&space;x_2\\&space;x_3&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;1\\&space;2\\&space;3&space;\end{bmatrix}" title="\begin{bmatrix} 1 &0 &0 \\ 0 &2 &0 \\ 0 &0 &3 \end{bmatrix}\begin{bmatrix} x_1\\ x_2\\ x_3 \end{bmatrix} = \begin{bmatrix} 1\\ 2\\ 3 \end{bmatrix}" /></a>

**Input:** 

*m* : INTEGER - number of rows and columns in the matrix *A*, and the length of *b* and *x*

*A* : REAL - square, diagonal matrix of size *m* x *m*

*b* : REAL - arbitrary vector of length *m*

**Output:** 

*x* : REAL - the solution vector of length *m*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: m, i
REAL*8, ALLOCATABLE :: A(:, :)
REAL*8, ALLOCATABLE :: x(:), b(:)

m = 3
ALLOCATE(A(1:m, 1:m), b(1:m), x(1:m))
A = RESHAPE((/1.D0, 0.D0, 0.D0, &
			& 0.D0, 4.D0, 0.D0, &
			& 0.D0, 0.D0, 7.D0/), (/m, m/), ORDER=(/2, 1/))
b = (/1.D0, 2.D0, 3.D0/)
CALL lss_diag(m, A, b, x)
WRITE(*,*) x
```

The outputs from the above code:

```fortran
   1.0000000000000000       0.50000000000000000       0.42857142857142855 
```

**Implementation/Code:** The code for lss_diag can be seen below.

```fortran
SUBROUTINE lss_diag(m, A, b, x)
	IMPLICIT NONE
	
	! Takes as an input a matrix `A` of size `m` x `m` and a vector `b`
	! of length `m`. Out puts a solution vector `x` of length `m` that
	! solves the system of equations `A``x`=`b`. This is a special case
	! where `A` is a diagonal matrix.
	INTEGER, INTENT(IN) :: m
	REAL*8, INTENT(IN) :: A(1:m, 1:m), b(1:m)
	REAL*8, INTENT(OUT) :: x(1:m)
	
	! Initialize an increment variable.
	INTEGER :: i
	
	! Find the solution to the linear system of equations by simply
	! dividing the value in `b` by the corresponding diagonal element
	! in `A`.
	DO i = 1, m
		x(i) = b(i)/A(i, i)
	END DO
END SUBROUTINE

```



