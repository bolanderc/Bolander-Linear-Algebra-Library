# Homework 8 Task 9: Rayleigh Quotient vs. Inverse Iteration Computational Time

*Task: Compare the Inverse Iteration Algorithm and the Rayleigh Quotient Algorithm in terms of the amount of time needed to compute an eigenvalue. Tabulate your results for the two methods.*

This process can be executed in Fortran with the following code

```fortran
INTEGER :: n, maxiter, i, j
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, t1, t0, alpha

n = 16

tol = 10.D-16
maxiter = 10000
alpha = 3.D0
ALLOCATE(A(n, n), v0(n), v(n))
CALL spd_mat_gen(A, n)
v0 = 1.D0
WRITE(*,*) "Inverse Iteration"
CALL CPU_TIME(t0)
CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, 1, v, lam)
CALL CPU_TIME(t1)
WRITE(*,*) lam
WRITE(*,*) t1 - t0
WRITE(*,*) "Rayleigh Quotient"
v0 = 1.D0
CALL CPU_TIME(t0)
CALL rayleigh_quotient(A, n, v0, tol, maxiter, 1, v, lam)
CALL CPU_TIME(t1)
WRITE(*,*) lam
WRITE(*,*) t1 - t0

```

The computational times, number of iterations, and final convergence errors for the Rayleigh Quotient method are shown in the table below for random matrices.

|      Matrix Size| Computational Time (s) |   # of Iterations   |
| :--: | :--: | :--: |
|   4   | 2.4999999999999849E-005 | 8 |
|   8   |   6.1999999999999989E-005   |   7   |
|    16  | 1.3499999999999992E-004 | 6 |
|     32 | 6.0799999999999960E-004 | 5 |

And for the Inverse Iteration method...

| Matrix Size | Computational Time (s) | # of Iterations |
| :---------: | :----------------: | :-------------: |
|      4       | 1.3200000000000017E-004 | 27 |
|       8      |         1.5700000000000002E-004           |       35       |
|        16     | 1.8699999999999988E-004 | 32 |
|         32    | 5.8299999999999975E-004 | 26 |

It should be noted that these are for a fixed guess for alpha, so there is some variation depending on how close the guess is to a given eigenvalue. In general, though, the trend follows that the Rayleigh Quotient and the Inverse Iteration take roughly the same amount of time (at times an order of magnitude different, though!), but the Rayleigh Quotient converges in FAR FEWER iterations.