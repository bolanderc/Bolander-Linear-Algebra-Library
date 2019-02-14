# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           sym_mat_gen

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c sym_mat_gen.f90```

and can be added to a program using

```$ gfortran program.f90 sym_mat_gen.o ``` 

**Description/Purpose:** This routine takes a square matrix and turns it into a symmetric matrix with values ranging from 0 to 1. For example, with a 3x3 matrix, the following diagram shows how the values would be linked.

**|x 1 2|**
**|1 x 3|**
**|2 3 x|**

where the **x**'s represent a random number that is not linked to any other and **1** links to **1** across the diagonal, etc. It can be used to quickly generate a symmetric matrix to be tested with other linear algebra routines.

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix

**Output:** 

*mat* : REAL - a symmetric array of size (n, n) containing random numbers

**Usage/Example:**

This routine has two inputs, `r` and `c`, and can be implemented in a program as follows

```fortran
n = 4

ALLOCATE(mat(n, n))
CALL sym_mat_gen(n, mat)
DO i = 1, n
WRITE(*, *) mat(i, :)
END DO
```

The outputs from the above code:

```fortran
  0.75751776322702891       0.49647460323199966       0.66837958508287576     
  0.49647460323199966       0.94740982034736176       0.77361628542110361     
  0.66837958508287576       0.77361628542110361       0.53126456636111374 
```

**Implementation/Code:** The code for sym_mat_gen can be seen below.

```fortran
SUBROUTINE sym_mat_gen(n, mat)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(OUT) :: mat(n, n)
INTEGER :: i, j

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
END DO

END SUBROUTINE
```



