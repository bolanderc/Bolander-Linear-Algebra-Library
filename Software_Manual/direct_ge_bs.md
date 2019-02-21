# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           direct_ge_bs

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c direct_ge_bs.f90```

and can be added to a program using

```$ gfortran program.f90 direct_ge_bs.o ``` 

**Description/Purpose:** This routine uses the subroutines `mat_row_ech` and `backsub` to solve a square linear system of equations such as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

using Gaussian elimination to reduce the augmented coefficient matrix to row echelon form and then apply the backward substitution method to find an approximate solution for ***x***. The augmented coefficient matrix is defined by

<a href="https://www.codecogs.com/eqnedit.php?latex=\left[&space;\begin{array}{ccc|c}&space;a_{11}&space;&&space;a_{12}&space;&&space;a_{13}&space;&&space;b_{1}&space;\\&space;0&space;&&space;a_{22}&space;&&space;a_{23}&space;&&space;b_{2}&space;\\&space;0&space;&&space;0&space;&&space;a_{33}&space;&&space;b_{3}&space;\\&space;\end{array}&space;\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[&space;\begin{array}{ccc|c}&space;a_{11}&space;&&space;a_{12}&space;&&space;a_{13}&space;&&space;b_{1}&space;\\&space;0&space;&&space;a_{22}&space;&&space;a_{23}&space;&&space;b_{2}&space;\\&space;0&space;&&space;0&space;&&space;a_{33}&space;&&space;b_{3}&space;\\&space;\end{array}&space;\right]" title="\left[ \begin{array}{ccc|c} a_{11} & a_{12} & a_{13} & b_{1} \\ 0 & a_{22} & a_{23} & b_{2} \\ 0 & 0 & a_{33} & b_{3} \\ \end{array} \right]" /></a>

The Gaussian elimination process computes the action of the inverse of ***A*** on ***b***.

**Input:** 

*m* : INTEGER - number of rows in the augmented coefficient matrix *A* (corresponds to the length of the solution vector *x*)

*n* : INTEGER - number of columns in the augmented coefficient matrix *A* (includes the square matrix x as well as the )

*aug_A* : REAL - augmented coefficient matrix of size *m* x *n*

**Output:** 

*x* : REAL - the solution to the system of equations in *aug_A*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, i
REAL*8, ALLOCATABLE :: A(:, :), x(:)

n = 4
m = 3
ALLOCATE(A(1:m, 1:n), x(1:m))
A = RESHAPE((/2.D0, 3.D0, 3.D0, -3.D0, &
			& 1.D0, -3.D0, 5.D0, 8.D0, &
			& 4.D0, 4.D0, 12.D0, 4.D0/), (/m, n/), ORDER=(/2, 1/))
CALL direct_ge_bs(A, m, n, x)
WRITE(*,*) x
```

The outputs from the above code:

```fortran
  -1.7999999999999994       -1.1000000000000003        1.2999999999999998   
```

**Implementation/Code:** The code for direct_ge_bs can be seen below.

```fortran
SUBROUTINE direct_ge_bs(aug_A, m, n, x)
	IMPLICIT NONE
	
	! Takes as inputs an augmented coefficient matrix, `aug_A` of size
	! `m` x `n` and outputs the solution when using Gaussian Elimination
	! and backward substitution, `x` of length `n`.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: aug_A(1:m, 1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Executes the `mat_row_ech` subroutine to take `aug_A` to row
	! echelon form.
	CALL mat_row_ech(aug_A, m, n)
	
	! Executes the `backsub` subroutine to find the solution to the
	! system of equations in `aug_A`.
	CALL backsub(m, aug_A(1:m, 1:n-1), aug_A(1:m, n), x)
	
	
END SUBROUTINE
```



