# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           abs_err_vecl2

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c abs_err_vecl2.f90```

and can be added to a program using

```$ gfortran program.f90 abs_err_vecl2.o l2_vec_norm.f90 ``` 

**Description/Purpose:** This routine calculates the error between an arbitrary vector ***a*** and its approximation, <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\hat{a}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\hat{a}}" title="\mathbf{\hat{a}}" /></a>, using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm or Euclidean norm of difference between the two, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=||\mathbf{\hat{a}}&space;-&space;\mathbf{a}||_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||\mathbf{\hat{a}}&space;-&space;\mathbf{a}||_2" title="||\mathbf{\hat{a}} - \mathbf{a}||_2" /></a>



**Input:**  

*n* : INTEGER - the length of the vector, *a*

*a* : REAL - an arbitrary vector of length *n*

*approx* : REAL - an approximation to the vector *a*

**Output:** 

*error* : REAL - the absolute error between *a* and *approx* using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 3
ALLOCATE(a(1:n), approx(1:n))
a(1) = 1
a(2) = 100
a(3) = 9
approx(1) = 1.1
approx(2) = 99
approx(3) = 11
CALL abs_err_vec(a, approx, n, norm)
WRITE(*,*) norm 
```

The output from the above code:

```fortran
   2.2383029296251151
```

which is the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm of the vector ***a***.

**Implementation/Code:** The code for l2_vec_norm can be seen. Note that another subroutine called in this 

```fortran
SUBROUTINE abs_err_vec(a, approx, n, error)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), approx(1:n)
REAL*8, INTENT(OUT) :: error

REAL*8 :: diff(1:n)
INTEGER :: i
REAL*8 :: norm

DO i = 1, n
	diff(i) = approx(i) - a(i)
END DO

CALL l2_vec_norm(diff, n, error)

END SUBROUTINE
```
