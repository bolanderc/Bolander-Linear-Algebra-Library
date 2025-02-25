# Calculate the Inner (Dot) Product of Two Vectors

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          out_prod_vec

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c out_prod_vec.f90```

and can be added to a program using

```$ gfortran program.f90 out_prod_vec.o ``` 

**Description/Purpose:** This routine will calculate the outer product, ***C***, of two vectors, ***a*** and ***b***, using

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{a}\otimes{\mathbf{b}}&space;=&space;\vv{\mathbf{C}}&space;=&space;\begin{bmatrix}&space;u_1v_1&space;&u_1v_2&space;&.&space;&.&space;&.&space;&u_1v_n&space;\\&space;u_2v_1&space;&u_2v_2&space;&&space;&&space;&&space;&.&space;\\&space;.&space;&.&space;&.&space;&&space;&&space;&.&space;\\&space;.&space;&.&space;&&space;&&space;&.&space;&.&space;\\&space;.&space;&.&space;&&space;&&space;&&space;&.&space;\\&space;u_mv_1&space;&u_mv_2&space;&.&space;&.&space;&.&space;&u_mv_n&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{a}\otimes{\mathbf{b}}&space;=&space;\vv{\mathbf{C}}&space;=&space;\begin{bmatrix}&space;u_1v_1&space;&u_1v_2&space;&.&space;&.&space;&.&space;&u_1v_n&space;\\&space;u_2v_1&space;&u_2v_2&space;&&space;&&space;&&space;&.&space;\\&space;.&space;&.&space;&.&space;&&space;&&space;&.&space;\\&space;.&space;&.&space;&&space;&&space;&.&space;&.&space;\\&space;.&space;&.&space;&&space;&&space;&&space;&.&space;\\&space;u_mv_1&space;&u_mv_2&space;&.&space;&.&space;&.&space;&u_mv_n&space;\end{bmatrix}" title="\mathbf{a}\otimes{\mathbf{b}} = \vv{\mathbf{C}} = \begin{bmatrix} u_1v_1 &u_1v_2 &. &. &. &u_1v_n \\ u_2v_1 &u_2v_2 & & & &. \\ . &. &. & & &. \\ . &. & & &. &. \\ . &. & & & &. \\ u_mv_1 &u_mv_2 &. &. &. &u_mv_n \end{bmatrix}" /></a>

**Input:**  

*m* : INTEGER - the length of vector *a*

*n* : INTEGER - the length of vector *b*

*a* : REAL - an arbitrary vector of length *m*

*b* : REAL - an arbitrary vector of length n

**Output:** 

*C* : REAL - the outer product of *a* and *b* with size *m* x *n*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, m, i
REAL*8, ALLOCATABLE :: a(:), b(:)
REAL*8, ALLOCATABLE :: C(:, :)

m = 4
n = 3
ALLOCATE(a(1:m), b(1:n), C(1:m, 1:n))
a = (/-3.2D0, 4.5D0, 8.7D0, 2.6D0/)
b = (/9.2D0, 4.6D0, -7.1D0/)
CALL out_prod_vec(m, n, a, b, C)
DO i = 1, m
	WRITE(*,*) C(i, :)
END DO
```

The output from the above code:

```fortran
   -29.439999999999998       -14.719999999999999        22.719999999999999     
   41.399999999999999        20.699999999999999       -31.949999999999999     
   80.039999999999992        40.019999999999996       -61.769999999999989     
   23.919999999999998        11.959999999999999       -18.460000000000001 
```

which is the vector outer product of ***a*** and ***b***.

**Implementation/Code:** The code for vec_cross_prod3 can be seen below.

```fortran
SUBROUTINE out_prod_vec(m, n, a, b, C)
	IMPLICIT NONE
	
	! Takes as inputs two vectors, `a` and `b`, of length `m` and `n`
	! respectively. Outputs the outer product of `a` and `b`, a matrix
	! `C` of size `m` x `n`.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: a(1:m), b(1:n)
	REAL*8, INTENT(OUT) :: C(1:m, 1:n)
	
	! Initialize two increment variables for looping.
	INTEGER :: i, j
	
	! Loop over all elements of `C` and multiply the appropriate
	! elements of `a` and `b` together to calculate the outer product.
	DO i = 1, m
		DO j = 1, n
			C(i, j) = a(i)*b(j)
		END DO
	END DO
END SUBROUTINE
```



