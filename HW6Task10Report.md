# Homework 6 Task 10: Stability of QR Factorization

*Task: Complete an internet search for sites that discuss the stability of various algorithms used in computing the QR factorization of both rectangular and square matrices. Give a brief description of what you found and include citations for the pages you find.*



## Stability of Different QR Factorization Methods 

There are several methods available to decompose a matrix into its Q and R factors, including the Gram-Schmidt orthogonalization methods, the Householder transformations, and the Givens rotations. There are several sources that speak to the stability of these individual methods that I will post here. First, the Gram-Schmidt (and its modified form) orthogonalization method are both known to be unstable<sup>[1]</sup>. However, the reason that it is unstable is made very clear when we remember that the Gram-Schmidt procedure at its core is a sequence of multiplications of A from the right by upper triangular matrices<sup>[2]</sup>. Those matrices can very easily be ill-conditioned! Even worse, if A is rectangular, the algorithm is very sensitive to rounding error since the algorithm determines whether each new column is linearly independent from the previous columns. The modified Gram-Schmidt algorithm was an attempt to stabilize the classical Gram-Schmidt algorithm, though it is not entirely successful<sup>[3]</sup>.

The Householder transformation remedies some of these issues by instead multiplying A from the left with orthogonal Householder matrices, which have a condition number of 1<sup>[2]</sup>. The matrix A can also be decomposed by applying orthogonal rotations, which is the case in the Givens rotations method. Compared to the Householder method, Givens rotations require roughly twice as many elementary operations, so Householder is faster, but the Givens rotations are slightly more accurate. In the end, they are both more stable than the Gram-Schmidt orthogonalization method.


* [Source 1](<http://www.math.iit.edu/~fass/477577_Chapter_4.pdf>)
* [Source 2](<https://www-old.math.gatech.edu/academic/courses/core/math2601/Web-notes/3num.pdf>)
* [Source 3](<http://people.inf.ethz.ch/gander/papers/qrneu.pdf>)