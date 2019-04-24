# Homework 7 Task 6: Conjugate Gradient on Hilbert Matrices

*Task: Try out the conjugate gradient method from the previous task on Hilbert matrices of order 4, 8, 16, and 32. Describe the results you obtained.*

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
    CALL cg_method(A, n, b, tol, maxiter, printit, x0)
    CALL abs_err_vecl2(x, x0, n, error)
    WRITE(*,*) "l2 absolute error"
    WRITE(*,*) error
END DO
```

The results of running this code are shown here

```fortran
 Hilbert           4 x           4
   6.5512287803440038E-036           6 ! Convergence error and number of iterations at exit.
 l2 absolute error          ! Between x = 1 and the approximation found by steepest_descent.
   5.2560823755821785E-014
 Hilbert           8 x           8
   1.5753844737292138E-029          13
 l2 absolute error
   3.4664412397916406E-005
 Hilbert          16 x          16
   2.1131287064811690E-029          21
 l2 absolute error
   3.8808602380047758E-005
 Hilbert          32 x          32
   5.3857943014331820E-030          34
 l2 absolute error
   3.6494733688589646E-005

```

As can be seen, the solution obtained by the conjugate gradient method is not even close to the given solution of *x* (equal to *n* ones) for Hilbert matrices greater than or equal to 8 x 8. This is because the Hilbert matrices are terribly ill-conditioned (as the columns and rows become more and more dependent on one another as the size is increased). The conjugate gradient method, however, is much more efficient and accurate than the steepest descent method, with results given [here](./HW7Task4Report.md).