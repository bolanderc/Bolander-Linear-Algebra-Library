# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          vec_cross_prod3

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c vec_cross_prod3.f90```

and can be added to a program using

```$ gfortran program.f90 vec_cross_prod3.o ``` 

**Description/Purpose:** This routine will calculate the cross product of two vectors, ***a*** and ***b*** of length 3, using

<a href="https://www.codecogs.com/eqnedit.php?latex=a\times{b}&space;=&space;(a_2b_3&space;-&space;a_3b_2)\vv{i}&space;-&space;(a_1b_3&space;-&space;a_3b_1)\vv{j}&space;&plus;&space;(a_1b_2&space;-&space;b_1a_2)\vv{k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a\times{b}&space;=&space;(a_2b_3&space;-&space;a_3b_2)\vv{i}&space;-&space;(a_1b_3&space;-&space;a_3b_1)\vv{j}&space;&plus;&space;(a_1b_2&space;-&space;b_1a_2)\vv{k}" title="a\times{b} = (a_2b_3 - a_3b_2)\vv{i} - (a_1b_3 - a_3b_1)\vv{j} + (a_1b_2 - b_1a_2)\vv{k}" /></a>

**Input:**  

*a* : REAL - an arbitrary vector of length 3

*b* : REAL - an arbitrary vector of length 3

**Output:** 

*crossprod* : REAL - the cross product of *a* and *b* (a vector of length 3)

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
a = (/2.D0, 3.D0, 4.D0/)
b = (/5.D0, 6.D0, 7.D0/)
CALL vec_cross_prod3(a, b, crossprod)
WRITE(*,*) crossprod
```

The output from the above code:

```fortran
  -3.0000000000000000        6.0000000000000000       -3.0000000000000000 
```

which is the vector cross product of ***a*** and ***b***.

**Implementation/Code:** The code for vec_cross_prod3 can be seen below.

```fortran
SUBROUTINE vec_cross_prod3(a, b, crossprod)
	IMPLICIT NONE
	
	! I/O variable initiation
	REAL*8, INTENT(IN) :: a(1:3), b(1:3)
	REAL*8, INTENT(OUT) :: crossprod(1:3)
	
	! Calculate the cross product for two vectors of length 3.
	crossprod(1) = a(2)*b(3) - a(3)*b(2)
	crossprod(2) = -(a(1)*b(3) - a(3)*b(1))
	crossprod(3) = a(1)*b(2) - a(2)*b(1)
END SUBROUTINE
```



