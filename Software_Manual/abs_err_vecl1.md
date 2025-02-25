# Calculate the Absolute Error of Two Vectors Using the L1 Norm

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           abs_err_vecl1

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c abs_err_vecl1.f90```

and can be added to a program using

```$ gfortran program.f90 abs_err_vecl1.o l1_vec_norm.o ``` 

**Description/Purpose:** This routine calculates the error between an arbitrary vector ***a*** and its approximation, <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\hat{a}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\hat{a}}" title="\mathbf{\hat{a}}" /></a>, using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm of the difference, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=||\mathbf{\hat{a}}&space;-&space;\mathbf{a}||_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||\mathbf{\hat{a}}&space;-&space;\mathbf{a}||_1" title="||\mathbf{\hat{a}} - \mathbf{a}||_1" /></a>

**Input:**  

*n* : INTEGER - the length of the vector, *a*

*a* : REAL - an arbitrary vector of length *n*

*approx* : REAL - approximation of the vector *a*

**Output:** 

*error* : REAL - the absolute error between *a* and *approx* using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm.

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
CALL abs_err_vecl1(a, approx, n, norm)
WRITE(*,*) norm 
```

The output from the above code:

```fortran
   3.1000000000000001     
```

which is the absolute error between ***a*** and ***approx*** using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm.

**Implementation/Code:** The code for abs_err_vecl1 can be seen. Note that the l1_vec_norm subroutine is also called

```fortran
SUBROUTINE abs_err_vecl1(a, approx, n, error)
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

!~ Calculate the l1 vector norm using the difference between the vectors
!~ to find the error.
CALL l1_vec_norm(diff, n, error)

END SUBROUTINE
```

