# Calculate the Product of a Scalar and a Vector

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           s_mult_vec

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c s_mult_vec.f90```

and can be added to a program using

```$ gfortran program.f90 s_mult_vec.o ``` 

**Description/Purpose:** This routine multiply a scalar, *s*, into a vector, ***a***, and calculate the product.

**Input:**  

*n* : INTEGER - the length of the vector being multiplied by the scalar

*a* : REAL - an arbitrary vector of length *n*

*s* : REAL - an arbitrary scalar

**Output:** 

*a* : REAL - the product of *s* and *a*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, i
REAL*8 :: s, doo
REAL*8, ALLOCATABLE :: a(:)

n = 8
s = 1.2D0
doo = 1.0D0
ALLOCATE(a(1:n))
DO i = 1, n
	a(i) = doo
END DO
CALL s_mult_vec(s, a, n)
DO i = 1, n
	WRITE(*,*) a(i)
END DO
```

The output from the above code:

```fortran
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000     
   1.2000000000000000 
```

which is the vector ***a*** multiplied by *s*.

**Implementation/Code:** The code for s_mult_vec can be seen below.

```fortran
SUBROUTINE s_mult_vec(s, a, n)
	IMPLICIT NONE
	
	! Takes a vector a, its length, n, and a scalar, s as inputs. s is
	! multiplied into a and then a is returned.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: a(1:n)
	REAL*8, INTENT(IN) :: s
	
	! An increment variable to loop over the vector a
	INTEGER :: i

	! Multiply a vector, a, by a scalar, s, element-wise.
	DO i = 1, n
		a(i) = s*a(i)
	END DO
END SUBROUTINE
```



