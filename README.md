(C) 2012, Marius Posta (mariusposta@gmail.com)
Check LICENSE.txt for the legal blah-blah.

Uncapacitated Facility Location Problem solver
==============================================


I. Building

In order to build the code (by invoking 'make all'), you'll need to get the 
source for BTT from Antonio Frangioni (frangio@di.unipi.it) and replace the 
appropriate files in the BTT subdirectory.

- Why not include BTT?
- Because it's freely available only under a weird academic license.

- Why not use another non-differentiable optimizer instead of BTT?
- Because it's the best one out there and I don't have the expertise or the
  motivation to write a better one myself (seriously looked into it, though).

- What's this other stuff in the subdirectory?
- It's a C99 wrapper for BTT written by Paul-Virak Khuong, a friend and 
  colleague who had to kindness to share it with me (thanks Paul!). To put it
  mildly, BTT is not easily amenable to integration with other programs.

Execution performance has been found to be slightly sensitive to the
compiler used, and as of now we recommend Clang over GCC. This recommendation
may already be out of date as both compilers are now being constantly improved.


II. Running 

Basically, the executable runs the following main loop:

0. read instance size from stdin
1. (optional) read resolution parameters from stdin
2. read instance prices from stdin
3. solve the UFLP, printing progress information in stderr
4. print best known solution in stdout
5. read 'quit' from stdin or go back to 1.

This loop structure exists because the UFLP often occurs as a subproblem within
a larger optimization algorithm, in which the prices differ often only slightly
from one iteration to the next.

Resolution parameters may also be provided as command-line arguments.
A detailed description can be obtained by invoking 'solver --help'.

Example: to solve an instance 'cap_a' stored as './instances/beasley/cap74.bz2'
with a time limit of 10 seconds, invoke the following command at a bash prompt:
  echo `bzcat ./instances/beasley/cap74.bz2` quit | ./solver -n cap_74 -t 10

  
