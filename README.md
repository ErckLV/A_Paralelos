**Parallel Algorithms- CC301**
============================

UNI - Facultad de Ciencias
2014 - 1 


**Taught by: Glen Rodriguez**

**Modified algorithms by : Jesus Lovon**


Aboute the code files:
-----------------------

- hola.c : Basic communication test between processors. Recognize master worker

**Basic Integration**:
- ejemplo1.c : Basic example of integration using rectangles, calculating pi using a quarter of a circle
- trapecio.c: Same as ejemplo1.c but using trapezoidal rule for integration.
- simpson.c: Same as ejemplo1.c but using Simpson's rule
- ejemplo2.c: Calculating PI using Monte Carlo Method, simple method with square.
- ejemplo3.c/ejemplo3-1.c: Calculating PI using Monte Carlo Method, better precision reducing area.
- ejemplo3-2.c: Calculating PI using Monte Carlo Method, better precision faster than ejemplo3.c using 2 areas.

**Linear Systems**:
- linear.c: Sequential program using Monte Carlo with random walk to solve a linear system
- linear_eq_p1.c: Parallel program doing the same as linear.c
- linear_eq_p2.c: Parallel Monte Carlo using a different estimator
- linear_eq_p3.c: Parallel Monte Carlo reusing calculated values
