SUBROUTINE cg_method(A, n, b, tol, maxiter, printit, x0)
	IMPLICIT NONE
	
	! Implements the conjugate gradient method for solving a linear
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

	! References
	! --------------------
	! Ascher, U., and C. Greif. "A First Course in Numerical Methods.
	! SIAM." Society for (2011).
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	! Sets up temporary variables that will be used in the algorithm.
	REAL*8 :: alpha, d0, d1, temp, bdel, tollim
	REAL*8 :: r0(n), r1(n), p0(n), p1(n), x1(n), s(n)
	INTEGER :: k
	
	! Initialize variables
	alpha = 0.D0
	d0 = 0.D0
	d1 = 0.D0
	temp = 0.D0
	bdel = 0.D0
	tollim = 0.D0
	r0 = 0.D0
	r1 = 0.D0
	p0 = 0.D0
	p1 = 0.D0
	x1 = 0.D0
	s = 0.D0
	k = 0
	
	! Calculate the initial residual using the initial guess `x0`
	CALL mat_prod(A, x0, n, n, 1, r0)
	r0 = b - r0
	
	! Set initial search direction equal to the residual.	
	p0 = r0
	
	! Calculate a limit in the convergence error according to the
	! algorithm outlined in Ascher.
	CALL vec_dot_prod(b, b, n, bdel)
	tollim = tol*tol*bdel
	
	! Calculates the initial "error" using the residual
	CALL vec_dot_prod(r0, r0, n, d0)
	
	! Iterates using the steepest descent algorithm until the "error" is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (d0 > tollim .AND. k < maxiter)
	
		! Calculate the product of `A` and the search direction `p0`.
		CALL mat_prod(A, p0, n, n, 1, s)
		
		! Calculate the step size, `alpha`.
		CALL vec_dot_prod(p0, s, n, temp)
		alpha = d0/temp
		
		! Increment `x0` and `r0` using the step size and the direction
		! `p0` that will drive the solution to a minimum error.
		x1 = x0 + alpha*p0
		r1 = r0 - alpha*s
		
		! Calculate the "error" `d1` (for delta).
		CALL vec_dot_prod(r1, r1, n, d1)
		
		! Update the search direction according to the current residual
		! and the previous search direction.
		p1 = r1 + (d1/d0)*p0
		
		! Update variables and increment the counter.
		x0 = x1
		r0 = r1
		p0 = p1
		d0 = d1
		k = k + 1
	END DO
	
	! If the `printit` argument is provided as 1, then print the final
	! convergence error and number of iterations.
	IF (printit == 1) THEN
		WRITE(*,*) d0, k
	END IF
	
END SUBROUTINE
	
