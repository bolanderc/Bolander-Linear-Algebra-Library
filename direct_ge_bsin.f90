SUBROUTINE direct_ge_bsin(aug_A, m, n, x)
	IMPLICIT NONE
	
	! Takes as inputs an augmented coefficient matrix, `aug_A` of size
	! `m` x `n` and outputs the solution when using Gaussian Elimination
	! and backward substitution, `x` of length `n`.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(INOUT) :: aug_A(1:m, 1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	
	! Initialize decrement variables and a variable to compute the sum
	! of previous solutions integrated into the algorithm for the 
	! backward substitution method.
	INTEGER :: k_b, j_b
	REAL*8 :: backsum
	
	! Initialize increment variables i, j, and k as well as a factor
	! variable to be used in the row echelon algorithm.
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	
	! Executes the `mat_row_ech` subroutine inline to take `aug_A` to
	! row echelon form.
	
	! Loop across all columns except for the last (never need to touch
	! that entry to cancel out what is below it). This targers the pivot
	! elements
	DO k = 1, n - 1
	
		! Loop through all of the rows except for the first one to make
		! them zeros using the algorithm. Makes all entries zero beneath
		! the pivot element.
		DO i = k + 1, m
		
			! Calculate the factor to reduce ith row value to zero
			factor = aug_A(i, k)/aug_A(k, k)
			
			! Loop through all columns to subtract the factor multiplied
			! by the previous row entry.
			DO j = k , n
				aug_A(i, j) = aug_A(i, j) - factor*aug_A(k, j)
			END DO
		END DO
	END DO
	
	! Executes the `backsub` subroutine inline to find the solution to
	! the system of equations in `aug_A`.
	
	! Calculate the last value in the solution vector `x`.
	x(m) = aug_A(m, n)/aug_A(m, n - 1)
	
	! Loop through the remaining rows in `x` to calculate the solution
	! using the backward substitution algorithm.
	DO k_b = m-1, 1, -1
		backsum = 0.D0
		DO j_b = k_b + 1, n
			backsum = backsum + aug_A(k_b, j_b)*x(j_b)
		END DO
		x(k_b) = (aug_A(k_b, n) - backsum)/aug_A(k_b, k_b)
	END DO
	
	
END SUBROUTINE
