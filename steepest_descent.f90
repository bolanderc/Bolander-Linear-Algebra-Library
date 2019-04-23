SUBROUTINE steepest_descent(A, n, b, tol, maxiter, printit, x0)
	IMPLICIT NONE
	
	! Implements the steepest descent algorithm for solving a linear
	! system of equations.
	
	! Takes as inputs the coefficient matrix `A` of rank `n` as well
	! as the right-hand side vector `b` of length `n`. The argument
	! `tol` is given as the tolerance of convergence error and the
	! `maxiter` argument gives a limit to the number of iterations in
	! the algorithm. Finally, `printit` prints to the screen the
	! convergence error and iteration count if set to a value of 1.
	! The `x0` argument is passed in as the initial guess of the
	! solution to the system of equations and is passed out as the
	! final approximation at the end of the subroutine.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	! Sets up temporary variables that will be used in the algorithm.
	REAL*8 :: alpha, r0(n), r1(n), x1(n), d, s(n), temp
	INTEGER :: k
	
	! Initialize variables
	k = 0
	r0 = 0.D0
	r1 = 0.D0
	x1 = 0.D0
	d = 0.D0
	temp = 0.D0
	
	! Calculate the initial residual using the initial guess `x0`
	CALL mat_prod(A, x0, n, n, 1, r0)
	r0 = b - r0
	
	! Calculates the initial "error" using the residual
	CALL vec_dot_prod(r0, r0, n, d)
	
	! Iterates using the steepest descent algorithm until the "error" is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (d > tol .AND. k < maxiter)
	
		! Calculate the product of `A` and the residual.
		CALL mat_prod(A, r0, n, n, 1, s)
		
		! Calculate the step size, `alpha`.
		CALL vec_dot_prod(r0, s, n, temp)
		alpha = d/temp
		
		! Increment `x0` and `r0` using the step size and the direction
		! `r0` that will drive the solution to a minimum residual.
		x1 = x0 + alpha*r0
		r1 = r0 - alpha*s
		
		! Calculate the "error" `d` (for delta).
		CALL vec_dot_prod(r1, r1, n, d)
		
		! Update variables and increment the counter.
		x0 = x1
		r0 = r1
		k = k + 1
	END DO
	
	! If the `printit` argument is provided as 1, then print the final
	! convergence error and number of iterations.
	IF (printit == 1) THEN
		WRITE(*,*) d, k
	END IF
	
END SUBROUTINE
		
