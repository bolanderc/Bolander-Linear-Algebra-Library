# Calculate the Dot Product of Two Vectors

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          vec_dot_prod

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c vec_dot_prod.f90```

and can be added to a program using

```$ gfortran program.f90 vec_dot_prod.o ``` 

**Description/Purpose:** This routine will calculate the dot product of two vectors, ***a*** and ***b*** of the same size *n*, using

<a href="https://www.codecogs.com/eqnedit.php?latex=a\cdot{b}&space;=&space;\sum_{i=1}^na_ib_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a\cdot{b}&space;=&space;\sum_{i=1}^na_ib_i" title="a\cdot{b} = \sum_{i=1}^na_ib_i" /></a>

**Input:**  

*n* : INTEGER - the length of the vectors for which the dot product is being calculated

*a* : REAL - an arbitrary vector of length *n*

*b* : REAL - an arbitrary vector of length *n*

**Output:** 

*dotprod* : REAL - the dot product of *a* and *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 3
ALLOCATE(a(1:n), b(1:n))
a(1) = 1
a(2) = 3
a(3) = -5
b(1) = 4
b(2) = -2
b(3) = -1
CALL vec_dot_prod(a, b, n, norm)
WRITE(*,*) norm
```

The output from the above code:

```fortran
   3.0000000000000000  
```

which is the vector dot product of ***a*** and ***b***.

**Implementation/Code:** The code for vec_dot_prod can be seen below.

```fortran
SUBROUTINE vec_dot_prod(a, b, n, dotprod)
	IMPLICIT NONE
	
	! I/O variables
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: a(1:n), b(1:n)
	REAL*8, INTENT(OUT) :: dotprod
	
	! Initialize subroutine variables
	INTEGER :: i
	
	dotprod = 0.D0
	
	! Calculate the dot product of the vectors a and b by multiplying
	! the elements together and summing them up
	DO i = 1, n
		dotprod = dotprod + a(i)*b(i)
	END DO
END SUBROUTINE
```



