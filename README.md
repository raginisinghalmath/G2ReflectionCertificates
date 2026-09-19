# Validated reflection shooting for nearly parallel G2 structures

This repository contains the computer-assisted part of the existence proof
for two cohomogeneity-one nearly parallel G2 structures obtained by reflection:

- the `tau4` solution on the Berger space;
- the `tau23` solution in the exactly U(1)-invariant subfamily `nu = 0`.

The code has four rigorous parts:

1. outward-rounded interval arithmetic for the regularized equations;
2. a validated regular-singular Taylor start, including parameter derivatives;
3. multiprecision CAPD C1 integration of the solution and variational equation;
4. a strict interval Krawczyk inclusion.

The scripts report success only after the singular start, the complete flow,
the real G2 chamber, the full interval Jacobian and the strict Krawczyk
inclusion have all been verified.

## Requirements

- Julia 1.12;
- a C++17 compiler;
- GMP and MPFR;
- CAPD compiled with multiprecision support.

The Julia environment fixes `IntervalArithmetic` at version `1.0.10`.

One way to install CAPD is:

```bash
git clone https://github.com/CAPDGroup/CAPD.git
cd CAPD
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="$HOME/capd" -DCAPD_ENABLE_MULTIPRECISION=ON
cmake --build . -j
cmake --install .
```

The Julia build script searches for CAPD in the following order:

1. the executable named by `CAPD_CONFIG`;
2. `$CAPD_PREFIX/bin/capd-config`;
3. `capd-config` on `PATH`;
4. the headers and library under `$CAPD_PREFIX`.

If `CAPD_PREFIX` is not set, it defaults to `$HOME/capd`. On Apple silicon,
GMP and MPFR are searched for under `/opt/homebrew`; this can be changed by
setting `GMP_PREFIX`.

## Set up the Julia environment

Run all commands from the repository root:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the Julia unit tests:

```bash
julia --startup-file=no --project=. test/runtests.jl
```

## Locate the approximate zeros

The following two programs use floating-point integration and damped Newton
iteration to locate approximate zeros of the shooting maps:

```bash
julia --startup-file=no --project=. scripts/approximate_S4.jl
julia --startup-file=no --project=. scripts/approximate_S23.jl
```

Both programs initially treat `(a,c,nu,T)` as four independent variables.
For `S23`, the default seed has `nu=0.005`; Newton iteration drives this
parameter numerically to zero.  Thus the enhanced U(1) symmetry is discovered
by solving the full four-equation shooting problem rather than imposed at the
approximation stage.  Optional initial guesses may be supplied on the command
line as `(a,c,nu,T)` for either map; run either program with `--help` for the
precise syntax.

These are non-rigorous computations.  They explain how the centres of the
validated parameter boxes are found, but they do not establish the existence
of exact zeros.  The numerical programs use the same regular-singular Taylor
orders and start times as the certificates: order `24` at `delta=1/8` for
`S4`, and order `8` at `delta=1/1024` for `S23`.  Exact existence and
uniqueness are established by the CAPD and Krawczyk programs below.

More precisely, each program constructs a high-order initial value from the
regularized equation, integrates the ordinary differential equation with an
adaptive Dormand--Prince method, and applies damped Newton iteration to the
shooting residual.  Positivity of the three chamber quantities is checked
throughout every right-hand-side evaluation and reported at the final time.
The source shared by the two programs is `src/ApproximateShooting.jl`; it does
not call CAPD and should be regarded only as a reproducible way to locate the
candidate zeros subsequently certified by the interval programs.  Having
observed numerically that the `S23` zero has `nu=0`, the rigorous computation
then works on that exact invariant slice and certifies the reduced
three-equation map.

## Search for a two-ended solution on S7

The additional program

```bash
julia --startup-file=no --project=. scripts/approximate_S7_gluing.jl
```

searches for parameters satisfying

```text
eta_(a_eta,c_eta,nu_eta)(T1) = xi_(a_xi,b_xi,c_xi)(T2).
```

All eight components of the two endpoint states are matched.  The two copies
of `a` and `c` are independent shooting parameters.  The left-hand `eta`
solution uses the existing high-order `(3,-1)` regular-singular recurrence.
The right-hand `xi` solution uses the leading `(1,-1)` regular-singular
expansion from the existing script and is integrated in its outward local
coordinate with `q_s = -Phi(q)`.  This reversed sign is what is required when
the global right-hand piece is written as `xi(T1+T2-t)`.  To reduce the
leading start error, the code integrates from both `delta_xi` and
`delta_xi/2` and Richardson-extrapolates the two endpoints.

With no arguments, the program starts from the exact homogeneous squashed-S7
parameters

```text
(a_eta,c_eta,nu_eta,T1) =
    (sqrt(5)/25, -3sqrt(5)/5, sqrt(5)/50, pi/6),

(a_xi,b_xi,c_xi,T2) =
    (sqrt(5)/25, -2sqrt(5)/25, 19sqrt(5)/25, pi/6).
```

A different initial guess can be supplied as eight numbers:

```bash
julia --startup-file=no --project=. scripts/approximate_S7_gluing.jl \
    a_eta c_eta nu_eta T1 a_xi b_xi c_xi T2
```

The change of variables used by Newton's method enforces

```text
a_eta > 0,  c_eta < 0,
a_xi > 0,   c_xi > a_xi,   a_xi*c_xi-b_xi^2 > 0,
T1 > delta_eta,   T2 > delta_xi.
```

The script prints the singular values and numerical rank of the transformed
eight-by-eight matching Jacobian.  Full rank means that the computed match is
locally isolated when `lambda` is fixed.  A positive-dimensional family would
require a genuine Jacobian kernel (or an additional varying normalization),
and the script prints a candidate kernel basis when the numerical rank is less
than eight.  Any apparent small singular value should be checked at higher precision.
As with the other approximation scripts, this program finds candidates but
does not prove the existence of an exact match.

Compile and test the CAPD driver:

```bash
julia --startup-file=no --project=. scripts/test_capd_setup.jl
```

## Run the two certificates

For the `tau4` certificate, run:

```bash
julia --startup-file=no --project=. scripts/certify_tau4_capd.jl
```

This validates a zero in the box centred at

```text
(0.042329982765393195,
 -0.06469940434457604,
  0.006759263945977987,
  0.7656914313714904)
```

with radius `1e-10` in every coordinate.

For the U(1)-invariant `tau23` certificate, run:

```bash
julia --startup-file=no --project=. scripts/certify_tau23_u1_capd.jl
```

Here `nu = 0` is imposed exactly. The reduced variables `(a,c,T)` are centred
at

```text
(0.025880416406250467,
 -0.07764124921873236,
  0.8578017401911906)
```

with radii `(1e-7,3e-7,5e-7)`.

A successful run ends with

```text
SUCCESS: strict Krawczyk inclusion certified with CAPD C1 flow.
```

Each certificate also prints the Julia and package versions and the SHA256
hashes of all source files used in the computation.

## Repository contents

```text
Project.toml
Manifest.toml
deps/g2_capd_flow.cpp
scripts/certify_tau4_capd.jl
scripts/certify_tau23_u1_capd.jl
scripts/test_capd_setup.jl
scripts/approximate_S4.jl
scripts/approximate_S23.jl
scripts/approximate_S7_gluing.jl
src/CapdFlow.jl
src/ApproximateShooting.jl
src/ApproximateS7Gluing.jl
src/G2CapdCertificates.jl
src/G2Equations.jl
src/IntervalKrawczyk.jl
src/SingularStart.jl
src/ValidatedArithmetic.jl
test/runtests.jl
```

The generated CAPD binary is placed in `deps/build/` and is not part of the
repository.
