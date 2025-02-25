# Generate a Symmetric, Diagonally-Dominant Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           sym_dd_mat_gen

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c sym_dd_mat_gen.f90```

and can be added to a program using

```$ gfortran program.f90 sym_dd_mat_gen.o ``` 

**Description/Purpose:** This routine takes a square matrix and turns it into a symmetric matrix with values ranging from 0 to 1. The matrix is then made diagonally dominant by summing up the values in each row and adding them to the diagonal element. More information on these two processes can be found in the documentation for [sym_mat_gen](./sym_mat_gen.md) and [mat_dd](./mat_dd.md).

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix

**Output:** 

*mat* : REAL - a symmetric, diagonally dominant array of size (n, n) containing random numbers

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: mat(:, :)

n = 3
ALLOCATE(mat(1:n, 1:n))
CALL sym_dd_mat_gen(n, mat)
DO i = 1, n
	WRITE(*,*) mat(i, :)
END DO
```

The outputs from the above code:

```fortran
   1.0150633661154935       0.24714688633236725       0.38542446285838927     
  0.24714688633236725        1.3062044887584361       0.86589152354173504     
  0.38542446285838927       0.86589152354173504        1.9965163149834613 
```

**Implementation/Code:** The code for sym_mat_gen can be seen below.

```fortran
SUBROUTINE sym_dd_mat_gen(n, mat)
	IMPLICIT NONE
	
	! Takes as input the size of the square, symmetric matrix, `n` and
	! the matrix `mat` to be created.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(n, n)
	
	! Initializes increment variables `i` and `j` as well as a `row_sum`
	! variable to add to the diagonal elements to make them dominant
	INTEGER :: i, j
	REAL*8 :: row_sum
	
	!~ Calculates a symmetric, square matrix by looping over the matrix
	DO i = 1, n
		!~ goes along all rows, i, and columns where j >= i
		DO j = i, n
			!~ A random number is assigned to the entries along the diagonal and to 
			!~ the right on each row (for example, a 3x3 matrix would assign random
			!~ values to all elements with an x in the following diagram
			!~ |x x x|
			!~ |- x x|
			!~ |- - x|
			CALL RANDOM_NUMBER(mat(i, j))
			!~ The random number assigned above is then copied across the diagonal.
			!~ The following diagram gives an example of how the values would match
			!~ |x 1 2|
			!~ |1 x 3|
			!~ |2 3 x|
			!~ The diagonal values are unique.
			mat(j, i) = mat(i, j)
		END DO
		row_sum = 0.D0
		! Loops through the columns and adds up all of the elements in
		! the ith row.
		DO j = 1, n
			row_sum = row_sum + mat(i, j)
		END DO
		! Adds the value of the row summation to the diagonal element.
		mat(i, i) = row_sum
	END DO
END SUBROUTINE
```



