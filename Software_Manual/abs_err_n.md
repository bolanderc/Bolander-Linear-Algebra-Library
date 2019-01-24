# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           abs_err_n

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c abs_err_n.f90```

and can be added to a program using

```$ gfortran program.f90 abs_err_n.o ``` 

**Description/Purpose:** This routine will compute the double precision absolute error between two numbers. This is given by 
$$
\begin{equation*}

e_{abs} = ||x - \bar{x}||

\end{equation*}
$$
where $x$ is the approximation and $\bar{x}$ is the actual value being approximated.

**Input:**  

$x$ : REAL - actual value, double precision

$y$ : REAL - approximation, double precision

**Output:** 

$e_{abs}$ : REAL - double precision absolute error between $x$ and $y$.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
x = 7.124645
y = 8.0
CALL abs_err_n(x, y, e_abs)
WRITE(*,*) e_abs
```

The output from the above code:

```fortran
  0.87535476684570312 
```

which is the absolute error between $x$ and $y$.

**Implementation/Code:** The code for abs_err_n can be seen [here](../abs_err_n.f90).

