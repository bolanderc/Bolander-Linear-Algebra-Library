# Calculate the L-Inf Norm of a Square Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          mat_infnorm

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c mat_infnorm.f90```

and can be added to a program using

```$ gfortran program.f90 mat_norm.o l1_vec_norm.o ``` 

**Description/Purpose:** This routine calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of an arbitrary square matrix, ***A***. This is done by finding the maximum row <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=||A||_\infty&space;=&space;\max_{1\leq{i}\leq{m}}\sum_{j=1}^n|a_{i\,j}|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||A||_\infty&space;=&space;\max_{1\leq{i}\leq{m}}\sum_{j=1}^n|a_{i\,j}|" title="||A||_\infty = \max_{1\leq{i}\leq{m}}\sum_{j=1}^n|a_{i\,j}|" /></a>  

**Input:**  

*n* : INTEGER - the number of rows and columns in the square matrix, *A*

*A* : REAL - a square matrix of size *n* x *n*

**Output:** 

*norm* : REAL - the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of the matrix, *A*.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 3
ALLOCATE(A(1:n, 1:n), Approx(1:n, 1:n))
A(1, 1) = 1D0
A(1, 2) = 2D0
A(1, 3) = 3D0
A(2, 1) = 4D0
A(2, 2) = 5D0
A(2, 3) = 6D0
A(3, 1) = 7D0
A(3, 2) = 8D0
A(3, 3) = 9D0
Approx(1, 1) = 0.44D0
Approx(1, 2) = 2.36D0
Approx(1, 3) = 3.04D0
Approx(2, 1) = 3.09D0
Approx(2, 2) = 5.87D0
Approx(2, 3) = 6.66D0
Approx(3, 1) = 7.36D0
Approx(3, 2) = 7.77D0
Approx(3, 3) = 9.07D0
CALL mat_1norm(A - Approx, n, norm)
WRITE(*,*) norm
```

The output from the above code:

```fortran
   2.4400000000000004     
```

which is the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of the matrix ***A***.

**Implementation/Code:** The code for mat_infnorm can be seen below. Note that the l1_vec_norm subroutine is also called

```fortran
SUBROUTINE mat_infnorm(A, n, norm)
	IMPLICIT NONE
	
	! I/O variables
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n)
	REAL*8, INTENT(OUT) :: norm
	
	! Initialize subroutine variables
	INTEGER :: i
	REAL*8 :: row_1norm
	
	norm = 0.D0
	row_1norm = 0.D0
	
	! Loop through each row of the matrix A and calculate the l1 norm
	! treating that row as a vector.
	DO i = 1, n
		CALL l1_vec_norm(A(i, :), n, row_1norm)
		! If a new maximum value is found from the row l1 norm, then
		! store it.
		IF (row_1norm > norm) THEN
			norm = row_1norm
		ENDIF
	END DO
END SUBROUTINE
```

