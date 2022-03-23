# One-Dimensional-Schr-dinger-Equation

Project written in C++ on basic simulations of eigenstates of the Hamiltonian in the One Dimensional Schrödinger Equation for common potential functions.

Please read the REPORT I wrote in the PDF file.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------

# Objectives and Process

I made this project in my Computational Physics class at
Fordham University. This was a research/simulation project that
I made to "investigate" the simplest cases of the one dimensional
time-independent Schrodinger's equation. 

  My goal was to implement a way of determining the eigenvalues
and eigenvectors of the Hamiltonian in those cases. To do so, I
based my work on what we saw in class on Givens Matrices.

  I then wrote a six pages report on what I found and what I was
able to accomplished in the time frame I was given. (This is the 
separate PDF one can find with my name on it)

----------------------------------------------------------------------------------------------------------------------------------------------------

# A couple of things about this project

- This project was written in C++ on a Linux platform and
  uses Gnuplot for all graph related matters.  

- The "Code.cpp" file is the main code in which the Schrodinger
  Equation is initialized in matrix and all the function 
  needed to diagonalize the Hamiltonian and graph the various
  eigenstate are called.

- The "matrixmath2.cpp" is the file that allows access to 
  the header file of the same name in which are called all 
  the functions related to matrix manipulations.

- The "Potential.cpp" is the file that allows access to 
  the header file of the same name in which the different 
  potential energy functions are defined.

- In the "TEST_eigen" folder is the code I used to verify 
  that the code was giving proper eigenvalues and vectors.

- In the "Results" folder are graphical results of the 
  eigenstates in each specific potential cases for both 50
  by 50 matrices and 100 by 100 matrices.

--------------------------------------------------------------------------------------------------------------------------------------

# A couple of things about the code

- If you want to compile it, please make sure you include 
  the two header files present in this directory (Matrix
  AND Potential). 

- In "Results," there are my results as previously mentioned.
  However, do not be alarmed if you do not see such a file 
  when you run the code. "DatFile" files should appear; they
  are the results.

- In the "Usage" lines there is a "n MAX." This is not a 
  limit on the code but on the number of graphs. Up to 100 
  graphs will be drawn, no matter the n selected. Nonetheless,
  if n > 100, the code, while taking longer, should run 
  properly but only show the 100 first graphs.
  I tried to use a for loop in c++, but gnuplot seemed not 
  to recognize the iterating factor. I then tried to create 
  a for loop inside the gnuplot command with " for [i=1:n] " 
  but this draw all the plots on one graph, making them 
  unreadable, hence the 100 lines.

- The Normalization Step of the code is not properly behaving
  for all potential functions at all energy states on the 
  boundaries. They tend to increase linearly instead of remaining
  at zero y-value.
  
  -------------------------------------------------------------------------------------------------------------------------------------------------
  -------------------------------------------------------------------------------------------------------------------------------------------------
  
  Fall 2021.
  
  R. Van Laer
