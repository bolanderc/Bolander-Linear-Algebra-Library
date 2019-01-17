# Homework 1 Task 10: Linear Algebra Libraries Report

*Task: Search the internet for sites that discuss linear algebra packages for solving linear algebra problems. Find a couple of sites that most closely line up with what you think we will be doing in class. Reference the sites in a brief discussion.*



## [LAPACK](http://www.icl.utk.edu/~mgates3/docs/lapack.html#linear)

LAPACK is a Fortran library that is used to solve linear algebra problems. Since I am new to coding in Fortran, this was interesting for me to see, since I have heard so much about its use. This to me seems like the ideal implementation of a linear algebra solver, since Fortran is a compiled language and is well suited to the kinds of heavy computational burdens that linear algebra requires. I was also interested to learn that LAPACK has some select functionality with the C language. Some common routines available in LAPACK that I am confident we will look into in this class include:

- Matrix norm routines
- LU linear systems solver
- Linear least squares
- Eigenvalue solver
- Basic Matrix manipulation (add, subtract, etc.)

## [SciPy](https://docs.scipy.org/doc/scipy/reference/linalg.html)

SciPy is the linear algebra library that is used in Python. I use this library pretty consistently in my research, and was interested in learning more about the background of the code as well as some of the functionality of which I may not be aware. I was very interested to learn that there is some overlapping between LAPACK and SciPy in that SciPy has some wrappers around LAPACK routines. The routines included in SciPy are routines that I expect we will touch on in some point in this class. Specifically, I noticed:

- norm() - Matrix or vector norm
- lstsq() - Compute least-squares solution to equation Ax=b
- eig() - Solve an ordinary or generalized eigenvalue problem of a square matrix

## [IMSL Numerical Libraries](https://www.roguewave.com/products-services/imsl-numerical-libraries)

The International Mathematics and Statistics (IMSL) Numerical Libraries are interesting because they actually span a large number of programming languages including functionality for C, Fortran, and Python.  It has largely the same type of capabilities as the above-mentioned libraries, but seeing its variability is interesting. Not only do these libraries do linear algebra-type routines, but they also do statistics and data mining as well.

