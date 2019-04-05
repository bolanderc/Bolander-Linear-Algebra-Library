SUBROUTINE qr_factor_gs(A, m, n, Q, R)
	IMPLICIT NONE
	
	! Takes as a input a coefficient matrix `A` of size `m` x `n` and
	! finds the QR factorization of that matrix using the classical
	! Gram-Schmidt orthogonalization algorithm. Outputs `Q` and `R`,
	! which are of size `m` x `n` and `n` x `n` respectively. Uses the
	! `vec_dot_prod` and `l2_vec_norm` subroutines.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: A(1:m, 1:n)
	REAL*8, INTENT(OUT) :: Q(1:m, 1:n), R(1:n, 1:n)
	
	INTEGER i, j
	
	! Iterates over each column to perform the QR factorization
	DO j = 1, n
		Q(:, j) = A(:, j)
		
		! Calculates the dot product of the columns to the left of the
		! current column to find R(i, j).
		DO i = 1, j - 1
			! This line is what makes this the "classical" method.
			CALL vec_dot_prod(A(:, j), Q(:, i), m, R(i, j))
			
			! Updates Q(j) using the previous columns and R(i, j)
			Q(:, j) = Q(:, j) - R(i, j)*Q(:, i)
		END DO
		
		! Computes the new Q column using R.
		CALL l2_vec_norm(Q(:, j), m, R(j, j))
		Q(:, j) = Q(:, j)/R(j, j)
	END DO
	
END SUBROUTINE
