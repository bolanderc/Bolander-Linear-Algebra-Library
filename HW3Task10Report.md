# Homework 3 Task 10: Frobenius Norm and Consistent Matrix Norms

*Task: Search the internet for sites that define and discuss the Frobenius matrix norm. Also, look for sites that define consistent matrix norms. Summarize in a paragraph or two what you found and as usual cite the sites.*



## [The Frobenius Norm](https://www.quora.com/What-is-the-significance-of-the-Frobenius-norm)

There were many sites that talked about how to define the Frobenius norm, but I was more curious about some of the merits of using the Frobenius norm, which is why I chose to talk about this Quora forum entry (though I recognize that forums may or may not be a great place to get information). Specifically, the answer by Garrett Thomas explains that one of the values of using the Frobenius norm is that it is useful for gradient-based methods, since the Frobenius norm is differentiable with respect to the individual entries of the matrix.

## [Consistent Norms](https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF-MAT3350/h07/undervisningsmateriale/chap13slides.pdf)

This slideshow goes through the definition for consistent matrix norms. The definition that it gives is ***Consistent Matrix Norms. A submultiplicative matrix norm which is defined for all m, n ∈ N, is said to be a consistent matrix norm.*** A submultiplicative matrix norm means that a matrix norm of a multiplication of two matrices should be bounded by the product of the norms of the two matrices, i.e.
<a href="https://www.codecogs.com/eqnedit.php?latex=||\vv{A}\vv{B}||&space;\leq&space;||\vv{A}||\,||\vv{B}||" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||\vv{A}\vv{B}||&space;\leq&space;||\vv{A}||\,||\vv{B}||" title="||\vv{A}\vv{B}|| \leq ||\vv{A}||\,||\vv{B}||" /></a>
To summarize, a consistent matrix norm is one which satisfies the above inequality for all sizes of matrices.