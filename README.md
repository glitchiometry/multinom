README

Requirements:
- basics.c/basics.h (from whw_clib)
- Patience

This project is a work in progress: very little of it has been tested yet. 
Multivariate polynomials are represented recursively (over the dimension) using arrays of void pointers (over degree)
To accommodate constant terms, the actual 'data' is accessed via a void pointer, which can be assigned to either
an integer, a double, an array_voidstar, an 'indeterminate', or NULL (which is or should be interpreted as identically zero,
but which could be parsed erroneously as an indeterminate constant in rarely used out-dated routines.)
Several basic operations are supported:
  - Arithmetic
  - Composition (of compatible multinomial functions)
Additionally, a 'multinomial evaluator' structure is included that can in principle facilitate calculations for certain
large/complex multinomials, by evaluating those multinomials in terms of simpler constituent multinomials when it is
possible to do so.  At present there does not yet exist a way to automatically determine which evaluator is optimal for
a given multinomial expression, but evaluators can at least be constructed concurrently with an associated multinomial.
The hope is that this feature will facilitate the discovery of efficient multivariate polynomial evaluation algorithms,
and eventually allow these algorithms to be studied phenomenologically/statistically.
There are several objectives for this project at this stage (e.g. testing, debugging, improving efficiency, adding
support for multinomials over arbitrary [commutative] rings, etc.), but my main priority is to implement basic algorithms
for computing certain algebraic characteristics of sets of multinomials (namely, given a list of multinomials, finding a
minimal set of independent terms; this would allow, for example, the possibility of characterizing the general form of
multinomials that can be evaluated with a given multinomial evaluator.)
