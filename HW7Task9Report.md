# Homework 7 Task 9: Conjugate Gradient vs. Jacobi Iteration on SPDs

*Task: Compare results on symmetric positive definite linear systems of equations using Jacobi versus the conjugate gradient methods. Disucss the results on systems with at least 500 equations and unknowns.*

This process can be executed in Fortran with the following code.

```fortran
INTEGER :: maxiter, printit, i, k
REAL*8, ALLOCATABLE :: A(:, :), b(:), xo(:), j(:), cg(:)
REAL*8 :: tol, error
INTEGER :: n(1:6)

n = (/500, 600, 700, 800, 900, 1000/)

DO k = 1, 6
    maxiter = 50000
    printit = 1
    IF (ALLOCATED(A)) THEN
    	DEALLOCATE(A, b, xo, j, cg)
    END IF
    ALLOCATE(A(n(k), n(k)), b(n(k)), xo(n(k)), j(n(k)), cg(n(k)))
    tol = 10.D-16
    error = 0.D0
    xo = 1.D0
    CALL spd_mat_gen(A, n(k))
    ! Makes the matrix diagonally dominant for Jacobi to work.
    DO i = 1, n(k)
    	A(i, i) = A(i, i)*DBLE(n(k))*5.D0
    END DO
    CALL mat_prod(A, xo, n(k), n(k), 1, b)
    cg = 0.D0
    CALL cg_method(A, n(k), b, tol, maxiter, printit, cg)
    CALL abs_err_vecl2(xo, cg, n(k), error)
    WRITE(*,*) error
    CALL jacobi_solve(A, n(k), b, tol, maxiter, j, printit)
    CALL abs_err_vecl2(xo, j, n(k), error)
    WRITE(*,*) error
END DO
```

The results of running this code are shown here

```fortran
  n =          500
 CG-method
   7.7471649058587555E-019          13
 l2-norm error   2.2366241246014454E-014
 Jacobi Iteration
   8.9509041826236192E-016          23
 l2-norm error   1.4081462427542381E-014
 

 n =          600
 CG-method
   1.0374833064379859E-017          12
 l2-norm error   2.8221292767733373E-014
 Jacobi Iteration
   6.4736570491389375E-016          20
 l2-norm error   2.4991364423678272E-014
 

 n =          700
 CG-method
   3.6285664726777449E-016          12
 l2-norm error   3.9403534090818864E-014
 Jacobi Iteration
   6.8438743594178853E-016          20
 l2-norm error   1.9951915554581547E-014
 

 n =          800
 CG-method
   1.4305809029658403E-017          12
 l2-norm error   4.8018525552153277E-014
 Jacobi Iteration
   5.9787339602818165E-016          20
 l2-norm error   3.8948104817008972E-014
 

 n =          900
 CG-method
   1.7011303202085758E-015          11
 l2-norm error   5.8707497585619614E-014
 Jacobi Iteration
   4.9650683064945462E-016          20
 l2-norm error   3.4865422370747881E-014
 

 n =         1000
 CG-method
   1.0110961605834339E-017          12
 l2-norm error   4.1889407683219731E-014
 Jacobi Iteration
   7.0216669371534024E-016          20
 l2-norm error   2.9025909921586275E-014


```

This task was interesting, because the Jacobi iteration and the CG-method really want two different types of matrices. Depending on how diagonally dominant the *A* matrix was, the Jacobi could be slightly faster than the CG-method. For example, with the diagonal multiplied by 50 instead of 5, the following results could be seen.

```fortran
 n =          500
 CG-method
   4.9146432091710194E-017          13
 l2-norm error   2.6688339058560935E-014
 Jacobi Iteration
   2.2204460492503131E-016          11
 l2-norm error   1.9482437089023217E-014
 

 n =          600
 CG-method
   1.4881260431208706E-014          12
 l2-norm error   3.4781712176990340E-014
 Jacobi Iteration
   9.9920072216264089E-016           9
 l2-norm error   1.8790012211091803E-014
 

 n =          700
 CG-method
   8.9623524540871118E-015          12
 l2-norm error   3.2994177111323028E-014
 Jacobi Iteration
   9.8052242617805960E-016           9
 l2-norm error   2.5088338340847902E-014
 

 n =          800
 CG-method
   1.1293505935689456E-015          12
 l2-norm error   3.5439946679252944E-014
 Jacobi Iteration
   9.8052242617805960E-016           9
 l2-norm error   2.7477963787545484E-014
 

 n =          900
 CG-method
   1.0659551770342064E-013          11
 l2-norm error   4.7665819416027150E-014
 Jacobi Iteration
   1.1102230246251565E-016          10
 l2-norm error   2.5973050138179879E-014
 

 n =         1000
 CG-method
   7.1384763991878202E-016          12
 l2-norm error   5.1575326664035426E-014
 Jacobi Iteration
   0.0000000000000000               10
 l2-norm error   4.0294672504184761E-014

```

In general, the CG-method is much more consistent, provided that a symmetric, positive definite matrix was provided.