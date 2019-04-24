# Homework 7 Task 4: Steepest Descent on Hilbert Matrices

*Task: Try out your steepest descent method on Hilbert matrices of size 4, 8, 16, 32. Explain your results.*

This process can be executed in Fortran with the following code.

```fortran
INTEGER :: n, i, j, k, maxiter, printit
REAL*8, ALLOCATABLE :: A(:, :), x(:), b(:), x0(:)
REAL*8 :: tol, error
INTEGER :: s(1:4)

maxiter = 50000
tol = 10.D-16
printit = 1
s = (/4, 8, 16, 32/)
DO k = 1, SIZE(s)
    IF(ALLOCATED(A)) THEN 
    DEALLOCATE(A, x, b, x0)
    END IF
    n = s(k)
    ALLOCATE(A(n, n), x(n), b(n), x0(n))
    WRITE(*,*) "Hilbert", n, "x", n
    DO i = 1, n
        DO j = 1, n
        	A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
        END DO
    END DO
    x = 1.D0
    CALL mat_prod(A, x, n, n, 1, b)
    x0 = 0.D0
    CALL steepest_descent(A, n, b, tol, maxiter, printit, x0)
    CALL abs_err_vecl2(x, x0, n, error)
    WRITE(*,*) "l2 absolute error"
    WRITE(*,*) error
END DO
```

The results of running this code are shown here

```fortran
 Hilbert           4 x           4
   9.9982091869048752E-016       34763 ! Convergence error and number of iterations at exit.
 l2 absolute error          ! Between x = 1 and the approximation found by steepest_descent.
   2.4328922231801894E-004
 Hilbert           8 x           8
   4.0705289325092022E-015       50000
 l2 absolute error
   3.3278565285813086E-003
 Hilbert          16 x          16
   2.5353592906581117E-014       50000
 l2 absolute error
   5.5038308709090744E-003
 Hilbert          32 x          32
   5.9635309136700095E-014       50000
 l2 absolute error
   8.5661352573237742E-003

```

As can be seen, the solution obtained by the steepest descent method is not even close to the given solution of *x* (equal to *n* ones). This is because the Hilbert matrices are terribly ill-conditioned (as the columns and rows become more and more dependent on one another as the size is increased). Even if the number of maximum iterations is increased, the l2 error does not improve by any significant amount.