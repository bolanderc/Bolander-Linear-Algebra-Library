SUBROUTINE K2_cond(A, n, v0, tol, maxiter, printit, cond)
	IMPLICIT NONE
	
	! Finds an approximation of the l2 condition number using the `power
	! method` and `inverse iteration` subroutines.
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, an initial guess for the eigenvector, `v0`,
	! a tolerance, `tol`, for exiting the iterative solver as well as a
	! maximum number of iterations, `maxiter`. Finally, a flag to print
	! the final number of iterations and convergence error, `printit` is
	! an input. The output of the subroutine is an approximation of the
	! l2 condition number. ***Note that `A` needs to be a symmetric
	! positive definite matrix.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), v0(n), tol
	REAL*8, INTENT(OUT) :: cond
	
	REAL*8 :: vhigh(n), vlow(n), lam_1, lam_n, alpha
	
	! Find the lowest eigenvalue using the inverse iteration method.
	alpha = 0.D0
	
	! Calculates the largest eigenvalue.
	CALL power_method(A, n, v0, tol, maxiter, printit, vhigh, lam_1)
	WRITE(*,*) lam_1
	
	! Calculates the smallest eigenvalue.
	CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, vlow, lam_n)
	WRITE(*,*) lam_n
	
	! Estimates the l2 condition number.
	cond = lam_1/lam_n
END SUBROUTINE
	
