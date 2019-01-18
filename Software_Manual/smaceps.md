# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           smaceps

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c smaceps.f90```

and can be added to a program using

```$ gfortran program.f90 smaceps.o ``` 

**Description/Purpose:** This routine will compute the single precision value for the machine epsilon. This is a routine for analyzing the behavior of any computer. This
usually will need to be run one time for each computer.

**Input:** There are no inputs needed.

**Output:** 

*mach_eps* : REAL - single precision machine epsilon value

*i* : INTEGER - the number of binary digits that represent machine epsilon

**Usage/Example:**

This routine has no inputs and can be implemented in a program as follows

```fortran
CALL smaceps(mach_eps, i)
WRITE(*,*) mach_eps, i
```

The outputs from the above code:

```fortran
   5.96046448E-08          24
```

The number of decimal digits that can be represented is roughly eight (E-8 on the
end of the first value) and the number of binary digits that represent machine epsilon are 24.

**Implementation/Code:** The code for smaceps can be seen [here](../smaceps.f90).

