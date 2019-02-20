SUBROUTINE mat_row_ech(A, m, n)
	IMPLICIT NONE
	
	! Takes as an argument a coefficient matrix `A`, of size `m`x`n`.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(INOUT) :: A(1:m, 1:n)
	
	! Initialize increment variables i, j, and k as well as a factor
	! variable to be used in the algorithm.
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! Loop across all columns except for the last (never need to touch
	! that entry to cancel out what is below it). This targers the pivot
	! elements
	DO k = 1, n - 1
	
		! Loop through all of the rows except for the first one to make
		! them zeros using the algorithm. Makes all entries zero beneath
		! the pivot element.
		DO i = k + 1, m
		
			! Calculate the factor to reduce ith row value to zero
			factor = A(i, k)/A(k, k)
			
			! Loop through all columns to subtract the factor multiplied
			! by the previous row entry.
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
		END DO
	END DO
END SUBROUTINE
