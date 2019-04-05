SUBROUTINE lu_factor(A, n)
	IMPLICIT NONE
	
	! Takes as inputs a square coefficient matrix, `A` of size
	! `n` x `n` and outputs the LU factorization of that matrix.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! Loops through all columns except for the last.
	DO k = 1, n - 1
		
		! Loops through the rows in the matrix `A`.
		DO i = k + 1, n
			
			! Calculates the factor to be used in both the Gaussian
			! Elimination algorithm used to find the U matrix and 
			! stored as the L matrix.
			factor = A(i, k)/A(k, k)
			
			! Performs the Gaussian Elimination and finds the U matrix
			! components.
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
			
			! Stores the factor used in the corresponding L matrix
			! component.
			A(i, k) = factor
		END DO
	END DO
END SUBROUTINE
