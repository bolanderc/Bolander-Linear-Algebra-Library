SUBROUTINE direct_ge_bs(aug_A, m, n, x)
	IMPLICIT NONE
	
	! Takes as inputs an augmented coefficient matrix, `aug_A` of size
	! `m` x `n` and outputs the solution when using Gaussian Elimination
	! and backward substitution, `x` of length `n`.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: aug_A(1:m, 1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Executes the `mat_row_ech` subroutine to take `aug_A` to row
	! echelon form.
	CALL mat_row_ech(aug_A, m, n)
	
	! Executes the `backsub` subroutine to find the solution to the
	! system of equations in `aug_A`.
	CALL backsub(m, aug_A(1:m, 1:n-1), aug_A(1:m, n), x)
	
	
END SUBROUTINE
