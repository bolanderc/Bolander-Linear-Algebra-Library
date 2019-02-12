# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual
This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           dmaceps

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c dmaceps.f90```

and can be added to a program using

```$ gfortran program.f90 dmaceps.o ``` 

**Description/Purpose:** This routine will compute the double precision value for the machine epsilon. This is a routine for analyzing the behavior of any computer. This
usually will need to be run one time for each computer.

**Input:** There are no inputs needed.

**Output:** 

*mach_eps* : REAL - double precision machine epsilon value

*i* : INTEGER - the number of binary digits that represent machine epsilon

**Usage/Example:**

This routine has no inputs and can be implemented in a program as follows

```fortran
CALL dmaceps(mach_eps, i)
WRITE(*,*) mach_eps, i
```

The outputs from the above code:

```fortran
   1.1102230246251565E-016          53
```

The number of decimal digits that can be represented is roughly sixteen (E-16 on the
end of the first value) and the number of binary digits that represent machine epsilon are 53.

**Implementation/Code:** The code for dmaceps can be seen below.

```fortran
SUBROUTINE dmaceps(mach_eps, i)
	IMPLICIT NONE
	
	REAL*8, INTENT(OUT) :: mach_eps
	INTEGER, INTENT(OUT) :: i
	REAL*8 :: one, approx
	
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
	
	WRITE(*,*) "Loop limit exceeded."
	
END SUBROUTINE
```



