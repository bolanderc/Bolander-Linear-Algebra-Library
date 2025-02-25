# Calculate the Sum of Two Vectors

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           vec_add

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c vec_add.f90```

and can be added to a program using

```$ gfortran program.f90 vec_add.o ``` 

**Description/Purpose:** This routine will add two vectors, ***a*** and ***b***, together and find the sum, a vector ***c***.

**Input:**  

*n* : INTEGER - the length of the vectors being added together

*a* : REAL - an arbitrary vector of length *n*

*b* : REAL - an arbitrary vector of length *n*

**Output:** 

*c* : REAL - the sum of *a* and *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 8
ALLOCATE(a(1:n), b(1:n), c(1:n))
DO i = 1, n
	a(i) = i
	b(i) = i + 1
END DO
CALL vec_add(a, b, c, n)
DO i = 1, n
	WRITE(*,*) c(i)
END DO
```

The output from the above code:

```fortran
   3.0000000000000000     
   5.0000000000000000     
   7.0000000000000000     
   9.0000000000000000     
   11.000000000000000     
   13.000000000000000     
   15.000000000000000     
   17.000000000000000 
```

which is the vector ***c***.

**Implementation/Code:** The code for vec_add can be seen below.

```fortran
SUBROUTINE vec_add(a, b, c, n)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), b(1:n)
REAL*8, INTENT(OUT) :: c(1:n)
INTEGER :: i

! Add two vectors together element-wise.
DO i = 1, n
	c(i) = a(i) + b(i)
END DO
END SUBROUTINE
```



