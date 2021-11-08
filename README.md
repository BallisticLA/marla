# Matlab Algorithms for Randomized Linear Algebra (MARLA)

MARLA is a Matlab library for prototyping algorithms in a future C++ library for randomized numerical linear algebra.
The library is meant to be "LAPACK-like" and organizes its functionality into high-level "drivers" and lower-level "computational routines".
Some of its functionality is currently untested -- this will be fixed early in the week of November 8.

MARLA has a companion Python package called [PARLA](https://github.com/BallisticLA/parla).
PARLA has an object-oriented design which is more flexible than the current implementations in PARLA.
We're in the process of standardizing a procedural API that's consistent across these two libraries.
The state of PARLA and MARLA's APIs and unit tests is summarized in [this Google Sheets spreadsheet](https://docs.google.com/spreadsheets/d/15vIS5wkaVB5lUoVQZqg7J_2qsK04ycVN47Mo2LuIKAo/edit?usp=sharing).
(It will be hard to read that spreadsheet without having the RandLAPACK design document on-hand.)