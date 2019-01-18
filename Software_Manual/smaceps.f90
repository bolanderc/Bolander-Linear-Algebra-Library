SUBROUTINE smaceps(mach_eps, i)
	IMPLICIT NONE
	
	REAL, INTENT(OUT) :: mach_eps
	INTEGER, INTENT(OUT) :: i
	REAL :: one, approx
	
	! Initialize variables. These are set an arbitrary number that we are
	! going to use to determine the machine precision of the computer.
	one = 1.0
	mach_eps = 1.0
	approx = one + mach_eps
	
	! Do loop to calculate machine precision. Divides the arbitrary
	! number chosen above by 2 at each iteration until the computer
	! cannot tell the difference between `approx` and `one`.
	DO i=1, 1000
		mach_eps = mach_eps / 2.
		approx = one + mach_eps
		
		! If the computer cannot tell the difference between `approx`
		! and `one`, the subroutine exits with `i` (the number of
		! iterations), and `mach_eps` (machine precision).
		IF(ABS(approx - one) == 0.0) RETURN
		
	END DO

END SUBROUTINE
	
