This folder contains the source codes related to the manuscript:

E. G. Birgin, J. L. Gardenghi, and A. Laurain, Asymptotic bounds on
the optimal radius when covering a set with minimum radius identical
balls, submitted.

All codes needed to reproduce the numerical experiments in the paper
are included. The third party codes Algencan 4.0.0 and GEOMPACK are
required.

ALGENCAN
========

To install Algencan 4.0.0, visit:

https://www.ime.usp.br/~egbirgin/sources/bmcomper/

(Note that Algencan 4.0.0 requires routines from BLAS and HSL as
well.)

It is assumed that Algencan 4.0.0 was successfully installed, that the
environment variable ALGENCAN points to the folder were Algencan 4.0.0
was installed, and that within folders
$ALGENCAN/sources/algencan/lib/, $ALGENCAN/sources/blas/lib/, and
$ALGENCAN/sources/hsl/lib/ there are files named libalgencan.a,
libblas.a, and libhsl.a, respectively.

GEOMETRY AND GEOMPACK
=====================

The two third party packages GEOMETRY and GEOMPACK are also
required. They correspond to files geometry.f90 and geompack2.f90,
taken from

https://people.sc.fsu.edu/~jburkardt/f_src/geometry/geometry.html

and

https://people.math.sc.edu/Burkardt/f_src/geompack2/geompack2.html

that are included in the current distribution. They are compiled
together with the other codes and no specific installation process is
required.

COVERING CODE
=============

This distribution includes:

(a) A file named covering.f90 that contains the main program with the
(i) regular polygon generator (from 3 up to 12 vertices) and (ii) the
three problems considered in the paper (Sketch of American map, Cesàro
Fractal, and Non-regular with holes). An additional file
modamerica.f90 borrowed from the distribution of Algencan 3.1.1 is
also included. The file modcesaro.f90 aim to support generation of
Cesàro Fractal problem.

(b) A makefile to compile the main program. Just type `make` to
compile. It was tested in a Linux operating system, and its
configuration may vary if another operating system is used.

(c) A script run-tests.sh that runs all case tests considered in the
paper.

(d) The Scharge random number generator is incldued in file drand.f90

(e) There is also a file named vorintpols.f90 that contains
subroutines to deal with Voronoi cells and convex polygons. It
includes the Sutherland-Hodgman algorithm for 2D clipping polygons and
an algorithm to compute the intersection of a convex polygon with a
circle.

HOW TO
======

To solve any of the problems considered in the manuscript, simple run
the script and check the output in the screen. To run the script to a
specific problem, type on terminal:

./run-tests.sh america|cesaro|stoyan|regpol [nvert]

It will run the specified problem with 10, 20, 30, ..., 100 balls and
a budget of 1000.