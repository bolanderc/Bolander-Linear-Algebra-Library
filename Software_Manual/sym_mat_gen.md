# Generate a Symmetric Matrix

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

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: mat(:, :)

n = 3

ALLOCATE(mat(n, n))
CALL sym_mat_gen(n, mat)
DO i = 1, n
WRITE(*, *) mat(i, :)
END DO
```

The outputs from the above code:

```fortran
  0.16192776212350513       0.87910223888966588       0.14078135957167615     
  0.87910223888966588       0.84431844571234449       0.35833527779431407     
  0.14078135957167615       0.35833527779431407       0.21247010336995686
```

**Implementation/Code:** The code for sym_mat_gen can be seen below.

```fortran
SUBROUTINE sym_mat_gen(n, mat)
	IMPLICIT NONE
	
	! Takes as an input the size of the matrix `n` and outputs a
	! symmetric, square matrix of that size.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(n, n)
	INTEGER :: i, j

	!~ Fills `mat` with random values
	CALL RANDOM_NUMBER(mat)

	!~ goes along all rows, i, and columns where j >= i
	!~ The random number assigned to the upper triangular elements is
	!~ then copied across the diagonal.
	!~ The following diagram gives an example of how the values would match
	!~ |x 1 2|
	!~ |1 x 3|
	!~ |2 3 x|
	!~ The diagonal values are unique.
	DO i = 1, n
		DO j = i, n
			mat(j, i) = mat(i, j)
		END DO
	END DO

END SUBROUTINE
```



