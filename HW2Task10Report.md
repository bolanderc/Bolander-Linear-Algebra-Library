# Homework 2 Task 10: Matrix Norms Report

*Task: Search the internet for sites that discuss matrix norms. Look for sites that explain induced matrix norms. Write a brief summary of what you find. Limit the discussion to no more than two or three paragraphs and include links to the sites you cite.*



## [Notes on Vector and Matrix Norms](http://www.cs.utexas.edu/users/flame/Notes/NotesOnNorms.pdf)

This set of course notes discusses the proofs associated with vector and matrix norms as well as touching on practical applications of norms in conditioning linear systems. It also covers the formulation of induced matrix norms. I found it very interesting to see the three conditions to defining a matrix norm met for each of the norms we discussed in class (since we didn't rigorously define them). I learned that norms are used to look at how error in the right hand of a system propagates into the solution to the system, and that the problem is that they ***do*** propagate.

## [Induced Norms](https://nptel.ac.in/courses/122104019/numerical-analysis/kadalbajoo/lec1/fnode3.html)

This website covers the basics of induced matrix norms, as well as giving a good summary of the L1, L2, and Linf induced matrix norms. It immediately presents an idea that I hadn't considered: that there is value in having consistency between how we calculate a vector and matrix norm. That is why matrix norms induced by vectors are valuable (rather than just by virtue of the fact that vector norms are much easier to calculate). The exception to this is the L2 induced matrix norm, which requires finding eigenvalues to calculate.
