SUBROUTINE vec_dot_prod(a, b, n, dotprod)
	IMPLICIT NONE
	
	! I/O variables
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: a(1:n), b(1:n)
	REAL*8, INTENT(OUT) :: dotprod
	
	! Initialize subroutine variables
	INTEGER :: i
	
	dotprod = 0.D0
	
	! Calculate the dot product of the vectors a and b by multiplying
	! the elements together and summing them up
	DO i = 1, n
		dotprod = dotprod + a(i)*b(i)
	END DO
END SUBROUTINE
