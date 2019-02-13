# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           mat_add

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c mat_add.f90```

and can be added to a program using

```$ gfortran program.f90 mat_add.o ``` 

**Description/Purpose:** This routine will add two matrices, ***A*** and ***B***, together and find the sum, a matrix ***C***.

**Input:**  

*n* : INTEGER - the number of rows in the matrices *A*, *B*, and *C*

*m* : INTEGER - the number of columns in the matrices *A*, *B*, and *C*

*A* : REAL - an arbitrary matrix of size *n* x *m*

*B* : REAL - an arbitrary matrix of size *n* x *m*

**Output:** 

*C* : REAL - the sum of *A* and *B*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, i
REAL*8, ALLOCATABLE :: A(:, :), B(:, :), C(:, :)

n = 3
m = 2
ALLOCATE(A(1:n, 1:m), B(1:n, 1:m), C(1:n, 1:m))
A = RESHAPE((/1.D0, 2.D0, &
			& 3.D0, 4.D0, &
			& 5.D0, 6.D0/), (/n, m/), ORDER=(/2, 1/))
B = RESHAPE((/1.D0, 2.D0, &
			& 3.D0, 4.D0, &
			& 5.D0, 6.D0/), (/n, m/), ORDER=(/2, 1/))
CALL mat_add(n, m, A, B, C)
DO i = 1, n
	WRITE(*,*) C(i, :)
END DO
```

The output from the above code:

```fortran
   2.0000000000000000        4.0000000000000000     
   6.0000000000000000        8.0000000000000000     
   10.000000000000000        12.000000000000000  
```

which is the matrix ***C***.

**Implementation/Code:** The code for mat_add can be seen below.

```fortran
SUBROUTINE mat_add(n, m, A, B, C)
	IMPLICIT NONE
	
	! Takes as inputs: two arrays, A and B of size n x m to be added
	! together. Outputs an array C that is the sum of A and B and has
	! size n x m. n represents the number of rows and m the number of
	! columns.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(IN) :: A(1:n, 1:m), B(1:n, 1:m)
	REAL*8, INTENT(OUT) :: C(1:n, 1:m)
	
	! Initializes increment variables i and j to loop over the matrices.
	INTEGER :: i, j
	
	! Loops over the matrix C and fills its elements with the sum of the
	! corresponding elements in A and B.
	DO i = 1, n
		DO j = 1, m
			C(i, j) = A(i, j) + B(i, j)
		END DO
	END DO
END SUBROUTINE
```



