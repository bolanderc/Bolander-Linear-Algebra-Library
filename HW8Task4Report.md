# Homework 7 Task 4: Steepest Descent on Hilbert Matrices

*Task: Test the code developed above (K2_cond) to compute the condition number of the Hilbert matrix as the size of the matrix is increased. Tabulate and/or graph the results from your work.*

This process can be executed in Fortran with the following code.

```fortran
INTEGER :: n, maxiter, i, printit, j, k
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, alpha, cond
INTEGER :: s(1:4)

tol = 10.D-18
maxiter = 10000
printit = 1
s = (/4, 8, 16, 32/)
DO k = 1, SIZE(s)
    IF(ALLOCATED(A)) THEN 
        DEALLOCATE(A, v0, v)
    END IF
    n = s(k)
    ALLOCATE(A(n, n), v(n), v0(n))
    WRITE(*,*) "Hilbert", n, "x", n
    DO i = 1, n
    	DO j = 1, n
    		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
    	END DO
    END DO
    v0 = 1.D0
    CALL K2_cond(A, n, v0, tol, maxiter, printit, cond)
    WRITE(*,*) cond
END DO
```

The results of running this code are shown here

|  Size of Hilbert Matrix    |   Condition Number   |
| :----: | :----: |
|   4x4   |   15513.738738930397   |
|   8x8   |  15257575817.032135    |
|   16x16   |   3.6910873798176121E+019   |
|32x32	| 1.1010765858137861E+018 	|


Hilbert matrices are ***terribly*** conditioned. This means that a small perturbation in the right hand side of a system of equations with these matrices will yield drastically different results.