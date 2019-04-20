SUBROUTINE gaussseidel_solve(A, n, b, tol, maxiter, x0)
	IMPLICIT NONE
	
	! Takes as an input the square, coefficient matrix `A` of size `n`
	! as well as the right-hand side vector `b`. In addition, a
	! tolerance argument `tol` is used to determine sufficient
	! convergence and a `maxiter` argument specifies the maximum number
	! of iterations before the solver exits. An initial guess for the
	! solution, `x0`, is given as an input. It stores the iterations of
	! the Gauss-Seidel approximation to the solution x and is output as
	! soon as the solver exits.
	INTEGER, INTENT(IN) :: n, maxiter
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	INTEGER :: i, j, k
	REAL*8 :: ax_sum, x1(n), error
	
	x1 = x0
	k = 0
	error = 10.D0*tol
	
	! Starts the iterative Gauss-Seidel solver.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Finds the next iteration of the vector x.
		DO i = 1, n
			k = k + 1
			ax_sum = b(i)
			
			! Sets the Gauss-Seidel routine apart from the Jacobi in
			! that it uses any updated values previously calculated in 
			! the iterations to give a more 'accurate' approximation.
			DO j = 1, i -1
				ax_sum = ax_sum - A(i, j)*x1(j)
			END DO
			DO j = i + 1, n
				ax_sum = ax_sum - A(i, j)*x0(j)
			END DO
			
			! The x vector approximation is updated, an error is
			! calculated, and the x0 vector is reset to the approximation
			! given in this iteration.
			x1(i) = ax_sum/A(i, i)
			CALL abs_err_vecl2(x0, x1, n, error)
			x1 = x0
		END DO
	END DO
END SUBROUTINE
