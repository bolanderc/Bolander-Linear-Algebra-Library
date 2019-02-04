# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           abs_err_vecl_inf

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c abs_err_vecl_inf.f90```

and can be added to a program using

```$ gfortran program.f90 abs_err_vecl_inf.o l_inf_vec_norm.o ``` 

**Description/Purpose:** This routine calculates the error between an arbitrary vector ***a*** and its approximation, <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\hat{a}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\hat{a}}" title="\mathbf{\hat{a}}" /></a>, using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of the difference, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=||\mathbf{\hat{a}}-\mathbf{a}||_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||\mathbf{\hat{a}}-\mathbf{a}||_\infty" title="||\mathbf{\hat{a}}-\mathbf{a}||_\infty" /></a>

**Input:**  

*n* : INTEGER - the length of the vector, *a*

*a* : REAL - an arbitrary vector of length *n*

*approx* : REAL - approximation of the vector *a*

**Output:** 

*error* : REAL - the absolute error between *a* and *approx* using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 3
ALLOCATE(a(1:n), approx(1:n))
a(1) = 1D0
a(2) = 100D0
a(3) = 9D0
approx(1) = 1.1D0
approx(2) = 99D0
approx(3) = 11D0
CALL abs_err_vecl_inf(a, approx, n, norm)
WRITE(*,*) norm 
```

The output from the above code:

```fortran
   2.0000000000000000      
```

which is the absolute error between ***a*** and ***approx*** using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm.

**Implementation/Code:** The code for abs_err_vecl_inf can be seen. Note that the l_inf_vec_norm subroutine is also called

```fortran
SUBROUTINE abs_err_vecl_inf(a, approx, n, error)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), approx(1:n)
REAL*8, INTENT(OUT) :: error

INTEGER :: i
REAL*8 :: diff(1:n)

!~ Find the difference between *a* and its approximation for each
!~ element.
DO i = 1, n
	diff(i) = approx(i) - a(i)
END DO

!~ Calculate the l_inf vector norm using the difference between the
!~ vectors to find the error.
CALL l_inf_vec_norm(diff, n, error)

END SUBROUTINE
```

