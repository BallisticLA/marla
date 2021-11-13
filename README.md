# Matlab Algorithms for Randomized Linear Algebra (MARLA)

MARLA is a Matlab library for prototyping algorithms in a future C++ library for randomized numerical linear algebra.
The library is meant to be "LAPACK-like" and organizes its functionality into high-level "drivers" and lower-level "computational routines".
All of its driver-level functionality has at least basic tests.

MARLA has a companion Python package called [PARLA](https://github.com/BallisticLA/parla).
PARLA has an object-oriented design which is more flexible than the current implementations in MARLA.
We're in the process of standardizing a procedural API that's consistent across these two libraries.
The state of MARLA and PARLA's APIs and unit tests is summarized in [this Google Sheets spreadsheet](https://docs.google.com/spreadsheets/d/15vIS5wkaVB5lUoVQZqg7J_2qsK04ycVN47Mo2LuIKAo/edit?usp=sharing).
(It will be hard to read that spreadsheet without having the RandLAPACK design document on-hand.)

## Notes on tests
* Tests are for basic correctness.
* We're using an object-oriented testing framework for eigendecomposition
  methods (evd1, evd2). These tests were written last. We might go
  back and use this as an organizing principle when cleaning up other tests.

## Important TODOs

* Clean up tests.
* Complete documentation (and revise existing documentation to follow
  a common style).

## Less important TODOs

 * Add tests for computational routines that aren't used in downstream
   drivers (refer to the spreadsheet for which routines those are).
 * Matlab doesn't have native support for abstract linear operators.
   So far this is only a minor limitation.
   The main way that PARLA uses abstract linear operators is in 
   a common interface for sketching operators that are only represented
   implicitly. In Matlab it would suffice to define a special class
   just for SRTTs. The special class would just need to overload
   "*" (matmul) and "'" (transpose).