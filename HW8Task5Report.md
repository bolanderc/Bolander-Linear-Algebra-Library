# Homework 8 Task 5: Eigenvalue Searching

*Task: Write a code that will compute approximations of the smallest and largest eigenvalues. Next subdivide the interval containing the smallest and largest eigenvalues and use the points in the subdivision to shift and locate other eigenvalues. Discuss your results.*

This process can be executed in Fortran with the following code using the [eigen_search](./Software_Manual/eigen_search.md) subroutine (additional examples are shown in the Software Manual entry) and a matrix with tightly clustered eigenvalues.

```fortran
INTEGER :: n, maxiter
REAL*8, ALLOCATABLE :: A(:, :), v0(:)
REAL*8 :: tol

tol = 10.D-18
maxiter = 10000
n = 3

ALLOCATE(A(n, n), v0(n))
A = RESHAPE((/ 1.D0, 0.D0, 0.D0, &
			 & 0.D0, 1.2D0, 0.D0, &
			 & 0.D0, 0.D0, 1.1D0/), (/n, n/), ORDER=(/2, 1/))
v0 = 1.D0
CALL eigen_search(A, n, v0, tol, maxiter, 5)
```

The results of running five searches are shown here

```fortran
 Minimum Eigenvalue:   1.0000000000000009     
 Maximum Eigenvalue:   1.1999999999999995     
 Search           1
    Alpha                     Eigenvalue
   1.1000000000000001                            NaN
 Search           2
    Alpha                     Eigenvalue
   1.0666666666666671        1.0000000000000000     
   1.1333333333333333        1.0000000000000002     
 Search           3
    Alpha                     Eigenvalue
   1.0400000000000007        1.0000000000000000     
   1.0800000000000005        1.0000000000000000     
   1.1200000000000003        1.0000000000000002     
   1.1600000000000001        1.0000000000000000     
 Search           4
    Alpha                     Eigenvalue
   1.0285714285714294        1.0000000000000000     
   1.0571428571428578        1.0000000000000000     
   1.0857142857142863        1.0000000000000000     
   1.1142857142857148        1.0000000000000004     
   1.1428571428571432        1.0000000000000000     
   1.1714285714285717        1.1999999999999906     
 Search           5
    Alpha                     Eigenvalue
   1.0222222222222230        1.0000000000000000     
   1.0444444444444452        1.0000000000000000     
   1.0666666666666673        1.0000000000000000     
   1.0888888888888895        1.0000000000000002     
   1.1111111111111116        1.0999999999999992     
   1.1333333333333337        1.0000000000000002     
   1.1555555555555559        1.0000000000000000     
   1.1777777777777780        1.2000000000000000
```

The inverse iteration method is shown to be inconsistent at identifying eigenvalues using shifting, though they are all eventually identified.