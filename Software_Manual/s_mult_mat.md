# Calculate the Product of a Scalar and a Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           s_mult_mat

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c s_mult_mat.f90```

and can be added to a program using

```$ gfortran program.f90 s_mult_mat.o ``` 

**Description/Purpose:** This routine multiply a scalar, *s*, into a matrix, ***A***, and calculate the product.

**Input:**  

*s* : REAL - an arbitrary scalar

*n* : INTEGER - the number of rows in the matrix *A*

*m* : INTEGER - the number of columns in the matrix *A*

*A* : REAL - an arbitrary matrix of size *n* x *m*

**Output:** 

*A* : REAL - the matrix representing the product of *s* and *A*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, i
REAL*8 :: s
REAL*8, ALLOCATABLE :: a(:, :)

n = 3
m = 2
s = 1.2D0
ALLOCATE(a(1:n, 1:m))
a = RESHAPE((/1.D0, 2.D0, &
			& 3.D0, 4.D0, &
			& 5.D0, 6.D0/), (/n, m/), ORDER=(/2, 1/))
CALL s_mult_mat(s, n, m, a)
DO i = 1, n
	WRITE(*,*) a(i, :)
END DO
```

The output from the above code:

```fortran
   1.2000000000000000        2.3999999999999999     
   3.5999999999999996        4.7999999999999998     
   6.0000000000000000        7.1999999999999993 
```

which is the matrix ***A*** multiplied by *s*.

**Implementation/Code:** The code for s_mult_mat can be seen below.

```fortran
SUBROUTINE s_mult_mat(s, n, m, A)
	IMPLICIT NONE
	
	! Takes as inputs a matrix A of size n x m and a scalar s. Outputs
	! the product of s and A as a modified version of A.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(IN) :: s
	REAL*8, INTENT(INOUT) :: A(1:n, 1:m)
	
	! Initializes two increment variables to loop over A.
	INTEGER :: i, j
	
	! Loops over every value in A and multiplies that value by s.
	DO i = 1, n
		DO j = 1, m
			A(i, j) = s*A(i, j)
		END DO
	END DO
END SUBROUTINE
```



