# Calculate the Absolute Error of Two Scalars

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           abs_err_n

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c abs_err_n.f90```

and can be added to a program using

```$ gfortran program.f90 abs_err_n.o ``` 

**Description/Purpose:** This routine will compute the double precision absolute error between two numbers. This is given by

 <a href="https://www.codecogs.com/eqnedit.php?latex=e_{abs}&space;=&space;||x&space;-&space;\bar{x}||" target="_blank"><img src="https://latex.codecogs.com/gif.latex?e_{abs}&space;=&space;||x&space;-&space;\bar{x}||" title="e_{abs} = ||x - \bar{x}||" /></a>

where *x* is the approximation and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\bar{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\bar{x}" title="\bar{x}" /></a> is the actual value being approximated.

**Input:**  

*x* : REAL - actual value, double precision

*y* : REAL - approximation, double precision

**Output:** 

*e_abs* : REAL - double precision absolute error between *x* and *y*.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
x = 7.124645
y = 8.0
CALL abs_err_n(x, y, e_abs)
WRITE(*,*) e_abs
```

The output from the above code:

```fortran
  0.87535476684570312 
```

which is the absolute error between *x* and *y*.

**Implementation/Code:** The code for abs_err_n can be seen below.

```fortran
SUBROUTINE abs_err_n(x, y, e_abs)
IMPLICIT NONE
REAL*8, INTENT(IN) :: x, y
REAL*8, INTENT(OUT) :: e_abs

! Calculates the absolute error between an approximation, y, and
! the actual value, x. These are single, double precision values.
e_abs = ABS(y - x)

END SUBROUTINE
```

