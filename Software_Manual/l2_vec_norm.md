# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           l2_vec_norm

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c l2_vec_norm.f90```

and can be added to a program using

```$ gfortran program.f90 l2_vec_norm.o ``` 

**Description/Purpose:** This routine calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm or Euclidean norm of an arbitrary vector, ***a***.

**Input:**  

*n* : INTEGER - the length of the vector, *a*

*a* : REAL - an arbitrary vector of length *n*

**Output:** 

*norm* : REAL - the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm of the vector, *a*.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 12
doo = 1.0D0
norm = 0.0D0
ALLOCATE(a(1:n))
DO i = 1, n
	a(i) = doo
END DO
CALL l2_vec_norm(a, n, norm)
WRITE(*,*) norm
```

The output from the above code:

```fortran
     3.4641016151377544 
```

which is the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm of the vector ***a***.

**Implementation/Code:** The code for l2_vec_norm can be seen.

```fortran
SUBROUTINE l2_vec_norm(a, n, norm)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(OUT) :: norm
REAL*8 :: summation
INTEGER :: i

summation = 0.0D0

! Sum up the squares of each value in the vector a.
DO i = 1, n
	summation = summation + (a(i)*a(i))
END DO

! Find the l2 norm of the vector by taking the square root of the sum
! of the squares
norm = SQRT(summation)
END SUBROUTINE
```



