# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           l_inf_vec_norm

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c l_inf_vec_norm.f90```

and can be added to a program using

```$ gfortran program.f90 l_inf_vec_norm.o ``` 

**Description/Purpose:** This routine calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of an arbitrary vector, ***a***.

**Input:**  

*n* : INTEGER - the length of the vector, *a*

*a* : REAL - an arbitrary vector of length *n*

**Output:** 

*norm* : REAL - the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of the vector, *a*.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 12
doo = 1.0D0
norm = 0.0D0
ALLOCATE(a(1:n))
DO i = 1, n
	a(i) = -doo
END DO
a(4) = 18.0D0
CALL l_inf_vec_norm(a, n, norm)
WRITE(*,*) norm
```

The output from the above code:

```fortran
   18.000000000000000  
```

which is the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of the vector ***a***.

**Implementation/Code:** The code for l_inf_vec_norm is found below.

```fortran
SUBROUTINE l_inf_vec_norm(a, n, norm)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(OUT) :: norm
INTEGER :: i

norm = 0.0D0

! The element in the vector with the maximum absolute value is found
! and represents the l_infinity norm.
DO i = 1, n
	IF (ABS(a(i)) > norm) THEN
		norm = ABS(a(i))
	ENDIF
END DO
END SUBROUTINE
```


