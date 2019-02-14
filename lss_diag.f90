SUBROUTINE lss_diag(m, A, b, x)
	IMPLICIT NONE
	
	! Takes as an input a matrix `A` of size `m` x `m` and a vector `b`
	! of length `m`. Out puts a solution vector `x` of length `m` that
	! solves the system of equations `A``x`=`b`. This is a special case
	! where `A` is a diagonal matrix.
	INTEGER, INTENT(IN) :: m
	REAL*8, INTENT(IN) :: A(1:m, 1:m), b(1:m)
	REAL*8, INTENT(OUT) :: x(1:m)
	
	! Initialize an increment variable.
	INTEGER :: i
	
	! Find the solution to the linear system of equations by simply
	! dividing the value in `b` by the corresponding diagonal element
	! in `A`.
	DO i = 1, m
		x(i) = b(i)/A(i, i)
	END DO
END SUBROUTINE
