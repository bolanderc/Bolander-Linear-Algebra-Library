SUBROUTINE ls_solvejacobi(A, b, m, n, x0, tol, maxiter, printit)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using a
	! modified Jacobi iteration method.
	
	! *NOTE* This algorithm was developed in conjunction with Jackson
	! T. Reid.
	INTEGER, INTENT(IN) :: m, n, maxiter, printit
	REAL*8, INTENT(IN) :: A(m, n), b(n), tol
	REAL*8, INTENT(OUT) :: x0(n)
	
	REAL*8 :: big_B(n, n), y(n), error, x1(n), sum_ax, alpha(n)
	INTEGER :: i, j, k
	
	x1 = 0.D0
	sum_ax = 0.D0
	error = 10.D0*tol
	k = 0
	alpha = 0.D0
	
	! Determine the alpha for each row by summing all of the row values
	! in `A`. Alpha is the parameter that ensures diagonal dominance.
	DO i = 1, n
		DO j = 1, n
			alpha(i) = alpha(i) + 2.D0*A(i, j)
		END DO
	END DO
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Solve using the modified Jacobi Iteration method.
	! First, add Ialpha to the B matrix
	DO i = 1, n
		big_B(i, i) = big_B(i, i) + alpha(i)
	END DO
	
	! Iteration loop for the algorithm. Iterates until the error is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Implements the modified Jacobi algorithm.
		DO i = 1, n
			! Starts with y = A^Tb.
			sum_ax = y(i)
			! Subtracts L_ATA*x_k
			DO j = 1, i - 1
				sum_ax = sum_ax - big_B(i, j)*x0(j)
			END DO
			! Adds alpha*I*x_k
			sum_ax = sum_ax + alpha(i)*x0(i)
			! Subtracts U_ATA*x_k
			DO j = i + 1, n
				sum_ax = sum_ax - big_B(i, j)*x0(j)
			END DO
			! Divides by the diagonal of A^T*A + alpha*I
			x1(i) = sum_ax/big_B(i, i)
		END DO
		
		! Increments the counter and calculates the absolute error
		! between x0 and x1 using the l2 norm. Then sets x0 = x1 for 
		! the next iteration.
		k = k + 1
		CALL abs_err_vecl2(x0, x1, n, error)
		x0 = x1
	END DO
	IF (printit == 1) THEN
		WRITE(*,*) error, k
	END IF
	
END SUBROUTINE
