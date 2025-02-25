# Calculate the Product of Two Matrices

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          mat_prod

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c mat_prod.f90```

and can be added to a program using

```$ gfortran program.f90 mat_prod.o ``` 

**Description/Purpose:** This routine will multiply to matrices, ***A*** and ***B***, of equal inner dimension together and find the product, a matrix ***C***.

**Input:**  

*n* : INTEGER - the number of rows in matrix *A*

*m* : INTEGER - the equal inner dimension of *A* and *B* (the number of columns in *A* and rows in *B*)

*p* : INTEGER - the number of columns in matrix *B*

*A* : REAL - an arbitrary matrix of size *n* x *m*

*B* : REAL - an arbitrary matrix of size *m* x *p*

**Output:** 

*C* : REAL - the matrix product of *A* and *B* of size *n* x *p*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, p
REAL*8, ALLOCATABLE :: A(:, :), B(:, :), C(:, :)

INTEGER :: i

n = 3
m = 2
p = 3
ALLOCATE(A(1:n, 1:m), B(1:m, 1:p), C(1:n, 1:p))
A = RESHAPE((/2.D0, 8.5D0, &
			& 4.D0, 2.1D0, &
			& 9.D0, 10.D0/), (/3, 2/), ORDER=(/2, 1/))
B = RESHAPE((/4.5D0, 8.6D0, 3.D0, &
			& 2.33D0, 4.22D0, 7.D0/), (/2, 3/), ORDER=(/2, 1/))
CALL mat_prod(A, B, n, m, p, C)
DO i = 1, n
	WRITE(*,*) C(i, :)
END DO
```

The output from the above code:

```fortran
   28.805000000000000        53.069999999999993        65.500000000000000     
   22.893000000000001        43.262000000000000        26.700000000000003     
   63.799999999999997        119.59999999999999        97.000000000000000     
```

which is the product of ***A*** and ***B***.

**Implementation/Code:** The code for mat_prod can be seen below.

```fortran
SUBROUTINE mat_prod(A, B, n, m, p, C)
	IMPLICIT NONE
	
	! Calculates the matrix product of A and B, which is C. A is of size
	! n x m and B is of size m x p, therefore C is of size n x p.
	INTEGER, INTENT(IN) :: n, m, p
	REAL*8, INTENT(IN) :: A(1:n, 1:m), B(1:m, 1:p)
	REAL*8, INTENT(OUT) :: C(1:n, 1:p)
	
	! Initialize i, j, and k as increments in do loops for the matrix
	! multiplication and also assign dotprodr_c to hold the value for
	! the dot product to be used in the calculation of the C element.
	INTEGER :: i, j, k
	REAL*8 :: dotprodr_c
	
	! Loop through the rows (i) and columns (j) of C and calculate the
	! value from the multiplication that belongs in that element.
	DO i = 1, n
		DO j = 1, p
			dotprodr_c = 0.D0
	
			! Calculate the dot product of the ith row of A and the jth
			! column of B and assign that value to C_i, j
			DO k = 1, m
				dotprodr_c = dotprodr_c + A(i, k)*B(k, j)
			END DO
			C(i, j) = dotprodr_c
		END DO
	END DO
END SUBROUTINE
```



