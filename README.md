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
* We're using an object-oriented testing framework for eigendecomposition
  methods (evd1, evd2). Depending on how that goes we might modify
  other tests to use a similar framework.
* Low-rank approximation tests sometimes fail randomly.
  This should only happen with tests that have an "_approx" suffix.
  We'll fix this once we take care of the TODO below about propagating
  random states through function calls. 

## Important TODOs

* Have each function take in a random state and use it for all
  random number generation. If it's not possible to write to the
  state so that it's changed on exit, then the function should
  be modified to return a random state in addition to its usual
  return values.

## Less important TODOs

 * Matlab doesn't have native support for abstract linear oeprators.
   So far this is only a minor limitation.
   The main way that PARLA uses abstract linear operators is in 
   a common interface for sketching operators that are only represented
   implicitly. In Matlab it would suffice to define a special class
   just for SRTTs. The special class would just need to overload
   "*" (matmul) and "'" (tranpose).