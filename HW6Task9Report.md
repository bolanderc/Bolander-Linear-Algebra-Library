# Homework 6 Task 9: Normal Equations vs. QR Factorization

*Task: Do an internet search for pages that discuss the difference between the solution of the least squares problem using the normal equations and the solution via QR factorization of the matrix. Make sure that you cite the sites you use and include a couple of paragraphs in your own words.*



## Differences Between the Methods

The solution of a least squares problem using the normal equations is relatively straight-forward. Assuming the directional derivatives of the residual are all zero, we can solve the normal equations, i.e. 
<a href="https://www.codecogs.com/eqnedit.php?latex=A^TAx&space;=&space;A^Tb" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A^TAx&space;=&space;A^Tb" title="A^TAx = A^Tb" /></a>

The issue with this method comes down to the stability of the algorithm. An ill-conditioned system will be even worse when solved with the normal equations<sup>[1]</sup>. If the condition number is large, then small relative changes to A can cause large changes to the span of A as well. The first source also mentions that the sensitivity of the least squares problem to perturbations in A can quickly be dominated by the condition number of A<sup>[1]</sup>.

On the other hand, the QR factorization is much more stable than the normal equations because it avoids forming A<sup>T</sup>A<sup>[2]</sup> . This comes at the cost of efficiency, as the QR factorization costs 2n<sup>2</sup>m - 2/3n<sup>3</sup> flops which is about twice the cost of the normal equations when m is large compared to n<sup>[3]</sup>. The stability of QR factorization makes it the more reasonable choice in many situations, since getting to a bad solution quickly is obviously not better than getting to a good solution more slowly. In addition, being able to find the QR factors of the A matrix can be useful in many other algorithms.


* [Source 1](<https://www.cs.cornell.edu/~bindel/class/cs3220-s12/notes/lec10.pdf>)
* [Source 2](<http://www.seas.ucla.edu/~vandenbe/133A/lectures/ls.pdf>)
* [Source 3](<https://sites.math.washington.edu/~morrow/498_13/demmelsvd.pdf>)