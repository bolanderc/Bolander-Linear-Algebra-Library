SUBROUTINE jacobi_solve(A, n, b, tol, maxiter, x0)
	IMPLICIT NONE
	
	! Takes as inputs the coefficient matrix, `A` of size `n` and the
	! right-hand side vector `b` of length `n`. Also as inputs are a
	! tolerance for exit of the algorithm, `tol` and a maximum number of
	! iterations `maxiter`. An initial guess, `x0` is input and refined
	! throughout the algorithm with each successive iteration. When the
	! algorithm exits, `x0` is output as the final approximation of x.
	INTEGER, INTENT(IN) :: n, maxiter
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	INTEGER :: i, j, k
	REAL*8 :: error, x1(n), sum_ax
	
	x1 = 0.D0
	sum_ax = 0.D0
	error = 10.D0*tol
	k = 0
	
	! Iteration loop for the algorithm. Iterates until the error is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Implements the Jacobi algorithm of x_{k+1} = x_k + D^{-1}*r_k.
		DO i = 1, n
			sum_ax = b(i)
			DO j = 1, i - 1
				sum_ax = sum_ax - A(i, j)*x0(j)
			END DO
			DO j = i + 1, n
				sum_ax = sum_ax - A(i, j)*x0(j)
			END DO
			x1(i) = sum_ax/A(i, i)
		END DO
		
		! Increments the counter and calculates the absolute error
		! between x0 and x1 using the l2 norm. Then sets x0 = x1 for 
		! the next iteration.
		k = k + 1
		CALL abs_err_vecl2(x0, x1, n, error)
		x0 = x1
	END DO
END SUBROUTINE
