# Homework 8 Task 8: Calculating the Condition Number

*Task: Look for internet sites that estimate the condition number of a matrix. Document the sites and as usual cite the web pages used in your explanation of the information found.*

## What is the Condition Number?

The condition number of a matrix reflects the sensitivity of the solution to changes in the right-hand side vector<sup>[1]</sup>. In other words, given a system of equations such that Ax = b, the condition number of A affects how sensitive x is to changes in b. The condition number is defined as 

<a href="https://www.codecogs.com/eqnedit.php?latex=\kappa(A)&space;=&space;||A||\cdot&space;||A^{-1}||" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\kappa(A)&space;=&space;||A||\cdot&space;||A^{-1}||" title="\kappa(A) = ||A||\cdot ||A^{-1}||" /></a><sup>[2]</sup>

Since estimating the inverse of a matrix can be very expensive (and difficult if the condition number is high!), this definition is rarely useful in computational linear algebra. Instead, other methods are used to estimate the condition number.

## Methods for Estimating the Condition Number

One way to find the condition number of a least squares problem is to perform a QR factorization. In this case, the condition number has an analytical solution relating to the R factor of the factorization<sup>[3]</sup>. Additional methods include finding the eigenvalues of a matrix and using them to help determine the condition number. One such approach uses a quasi Monte Carlo method to do so<sup>[4]</sup>. Ultimately, there are many ways to come up with an estimate of the condition number, though actually inverting the matrix is generally not an efficient way to do so.

[Source 1](<https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/>)

[Source 2](<http://www.cse.iitd.ernet.in/~dheerajb/CS210_lect07.pdf>)

[Source 3](<http://www.netlib.org/lapack/lawnspdf/lawn273.pdf>)

[Source 4](<http://www.imar.ro/journals/Mathematical_Reports/Pdfs/2013/3/4.pdf>)