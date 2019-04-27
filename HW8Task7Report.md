# Homework 8 Task 7: Rayleigh Quotient Condition Number

*Task: Use the Rayleigh Quotient algorthm to compute an approximation for the condition of the matrix. Test your code on a Hilbert matrix of reasonable size. Discuss the results from the work.*

This process can be executed in Fortran with the following code using the [rayleigh_cond](./Software_Manual/rayleigh_cond.md) subroutine (additional examples are shown in the Software Manual entry)

```fortran
INTEGER :: n, maxiter
REAL*8, ALLOCATABLE :: A(:, :)
REAL*8 :: tol, cond

n = 3

tol = 10.D-16
maxiter = 10000

ALLOCATE(A(n, n))
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL rayleigh_cond(A, n, tol, maxiter, cond)
WRITE(*,*) cond

```

The condition number of this rank 3 Hilbert matrix are shown below

```fortran
    524.05677758606339 
```

As expected, the condition number for this matrix is incredibly large. Now, for a 10 x 10 matrix the condition number is

```fortran
 Maximum Eigenvalue:   1.6959389969219503     
 Minimum Eigenvalue:   1.1115389158721985E-010
   15257576434.840221
```

which is to be expected as well. This is verifiable using [this paper](<https://www.ams.org/journals/mcom/1967-21-099/S0025-5718-1967-0223075-0/S0025-5718-1967-0223075-0.pdf>).

