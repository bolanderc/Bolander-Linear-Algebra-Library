# Homework 5 Task 8: QR factorization Using Gram-Schmidt on Hilbert Matrices

*Task: Try out your QR-factorization method from the previous task on the Hilbert matrices of sizes n=4,6,8,10. Test to see if the orthogonal matrix is really orthogonal by multiplying QtQ and comparing the result to the identity matrix. Explain the results you obtain.*

The output of running the Hilbert matrices can be seen [here](./HilbertQR.txt). As can be seen, the QtQ section of each of the successive matrices becomes less and less similar to the identity matrix as the size of the matrix increases. This essentially means that Q is becoming less and less orthogonal due to round off error. If this same process were to be implemented using the classical Gram-Schmidt implementation, this would likely be even worse. This is where the Gram-Schmidt algorithm becomes less viable, as it is less stable for these matrices with small values in them where round off error is more likely to be significant.
