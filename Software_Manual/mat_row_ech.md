# Find the Row Echelon Form of a Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           mat_row_ech

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c mat_row_ech.f90```

and can be added to a program using

```$ gfortran program.f90 mat_row_ech.o ``` 

**Description/Purpose:** This routine implements a method that performs elementary row operations on a matrix to take the matrix to row echelon form. Row echelon form is defined using two criteria:

* all nonzero rows (rows with at least one nonzero element) are above any rows of all zeroes (all zero rows, if any, belong at the bottom of the matrix), and
* the leading coefficient (the first nonzero number from the left, also called the pivot of a nonzero row is always strictly to the right of the leading coefficient of the row above it ([citation](https://en.wikipedia.org/wiki/Row_echelon_form)).

For example, with a linear system of equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

***A*** will be taken to an upper triangular form through the rows.

**Input:** 

*m* : INTEGER - number of rows in the matrix *A*

*n* : INTEGER - number of columns in the matrix *A*

*A* : REAL - arbitrary matrix of size *m* x *n*

**Output:** 

*A* : REAL - the transformed, row echelon form of the input matrix *A*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, i
REAL*8, ALLOCATABLE :: A(:, :)

n = 3
m = 3
ALLOCATE(A(1:m, 1:n))
A = RESHAPE((/2.D0, 3.D0, 3.D0, &
			& 1.D0, -3.D0, 5.D0, &
			& 4.D0, 4.D0, 12.D0/), (/m, n/), ORDER=(/2, 1/))
CALL mat_row_ech(A, m, n)
WRITE(*,*) "A:"
DO i = 1, m
	WRITE(*,*) A(i, :)
END DO
```

The outputs from the above code:

```fortran
 A:
   2.0000000000000000        3.0000000000000000        3.0000000000000000     
   0.0000000000000000       -4.5000000000000000        3.5000000000000000     
   0.0000000000000000        0.0000000000000000        4.4444444444444446  
```

**Implementation/Code:** The code for mat_row_ech can be seen below.

```fortran
SUBROUTINE mat_row_ech(A, m, n)
	IMPLICIT NONE
	
	! Takes as an argument a coefficient matrix `A`, of size `m`x`n`.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(INOUT) :: A(1:m, 1:n)
	
	! Initialize increment variables i, j, and k as well as a factor
	! variable to be used in the algorithm.
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! Loop across all columns except for the last (never need to touch
	! that entry to cancel out what is below it). This targers the pivot
	! elements
	DO k = 1, n - 1
	
		! Loop through all of the rows except for the first one to make
		! them zeros using the algorithm. Makes all entries zero beneath
		! the pivot element.
		DO i = k + 1, m
		
			! Calculate the factor to reduce ith row value to zero
			factor = A(i, k)/A(k, k)
			
			! Loop through all columns to subtract the factor multiplied
			! by the previous row entry.
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
		END DO
	END DO
END SUBROUTINE
```



