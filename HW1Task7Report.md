# Homework 1 Task 7: Parallel Processing Report

Running the compiled [hello](hello.f) file creates the output

	christian@Joberta:~/Python Projects/math5610$ ./hello
	hello world from thread           0
	hello world from thread           6
	hello world from thread           2
	hello world from thread           5
	hello world from thread           4
	hello world from thread           3
	hello world from thread           1
	hello world from thread           7
	There are           8  threads!

I assumed that the processors would finish in order, but was surprised that they didn't, and now I understand that since they are working in parallel they just print in the order that they finish.
