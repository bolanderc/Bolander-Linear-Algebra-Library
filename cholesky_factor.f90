SUBROUTINE cholesky_factor(A, n)
	IMPLICIT NONE
	
	! Takes as inputs a symmetric, positive definite matrix `A` of size
	! `n` and outputs a modified version of `A` with the Cholesky
	! factorization stored in the lower triangular and diagonal parts of
	! the matrix and the original elements in the upper triangular part.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! The Cholesky factorization method first loops through all of the
	! diagonal elements and takes their square root...
	DO k = 1, n - 1
		A(k, k) = SQRT(A(k, k))
		
		! Then loops through the lower triangular components to compute
		! the Cholesky decomposition.
		DO i = k + 1, n ! rows
			A(i, k) = A(i, k)/A(k, k)
		END DO
		DO j = k + 1, n ! rows
			DO i = j, n ! columns
				A(i, j) = A(i, j) - A(i, k)*A(j, k)
			END DO
		END DO
	END DO
	! The last diagonal value is simply factored to its square root.
	A(n, n) = SQRT(A(n, n))
	
END SUBROUTINE
