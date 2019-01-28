# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           s_mult_vec

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c s_mult_vec.f90```

and can be added to a program using

```$ gfortran program.f90 s_mult_vec.o ``` 

**Description/Purpose:** This routine multiply a scalar, *s*, into a vector, ***a***, and calculate the product, ***c***.

**Input:**  

*n* : INTEGER - the length of the vector being multiplied by the scalar

*a* : REAL - an arbitrary vector of length *n*

*s* : REAL - an arbitrary scalar

**Output:** 

*c* : REAL - the product of *s* and *a*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 8
s = 1.2D0
doo = 1.0
ALLOCATE(a(1:n), c(1:n))
DO i = 1, n
	a(i) = doo
END DO
CALL s_mult_vec(s, a, c, n)
DO i = 1, n
	WRITE(*,*) c(i)
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

which is the vector ***c***.

**Implementation/Code:** The code for s_mult_vec can be seen [here](../s_mult_vec.f90).

