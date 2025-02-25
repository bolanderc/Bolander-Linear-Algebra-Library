# Create a Diagonally-Dominant Square Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           mat_dd

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c mat_dd.f90```

and can be added to a program using

```$ gfortran program.f90 mat_dd.o ``` 

**Description/Purpose:** This routine takes a square matrix, ***A***, and fills each element of the matrix with a random number from 0 to 1 not including one. This matrix is then modified so that it is diagonally dominant, which means that

<a href="https://www.codecogs.com/eqnedit.php?latex=|a_{i\,i}|&space;\geq&space;\sum_{\substack{j=1&space;\\&space;j&space;\neq&space;i}}^n|a_{i\,j}|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?|a_{i\,i}|&space;\geq&space;\sum_{\substack{j=1&space;\\&space;j&space;\neq&space;i}}^n|a_{i\,j}|" title="|a_{i\,i}| \geq \sum_{\substack{j=1 \\ j \neq i}}^n|a_{i\,j}|" /></a>

It can be used to quickly generate a matrix to be tested with other linear algebra routines.

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix

**Output:** 

*mat* : REAL - a diagonally dominant matrix of size (n, n) containing random numbers

**Usage/Example:**

This routine has an input *n* and can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8, ALLOCATABLE :: A(:, :)

n = 3
ALLOCATE(A(1:n, 1:n))
CALL mat_dd(n, A)
DO i = 1, n
	WRITE(*,*) A(i, :)
END DO
```

The outputs from the above code:

```fortran
   1.6853128504659201       0.65315654864851747       0.13964096814871718     
   3.8754568170421222E-002  0.97335793596721065       0.23884389800574723     
  0.34318290516739380       0.75534218933621855        2.0417728343013199   
```

**Implementation/Code:** The code for mat_dd can be seen below.

```fortran
SUBROUTINE mat_dd(n, mat)
	IMPLICIT NONE
	
	! Takes in a value for the number of rows and columns in the
	! diagonally dominant matrix and outputs the matrix.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(1:n, 1:n)
	
	! Intialize do loop incrementers and a row summation variable to
	! help enforce diagonal dominance.
	INTEGER :: i, j
	REAL*8 :: row_sum
	
	! Fills the matrix with random numbers in all elements from 0 to 1
	! non-inclusive.
	CALL RANDOM_NUMBER(mat)
	
	! Loops through all rows and calculates a value for the summation
	! of all elements in that row.
	DO i = 1, n
		row_sum = 0.D0
		DO j = 1, n
			row_sum = row_sum + mat(i, j)
		END DO
		! Adds the value of the row summation to the diagonal element.
		mat(i, i) = row_sum
	END DO
END SUBROUTINE
```



