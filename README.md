# Validated reflection shooting for nearly parallel G2 structures

This repository contains the computer-assisted part of the existence proof
for two cohomogeneity-one nearly parallel G2 structures obtained by reflection:

- the `tau4` solution on the Berger space;
- the `tau23` solution in the exactly U(1)-invariant subfamily `nu = 0`.

Both certificates use `lambda = 6/sqrt(5)`. This release also contains the
floating-point searches listed below. The general-`k` experiments for
`k = 2,...,8` discussed in Remark 5.10 of the manuscript are not included.

The code has four rigorous parts:

1. outward-rounded interval arithmetic for the regularized equations;
2. a validated regular-singular Taylor start, including parameter derivatives;
3. multiprecision CAPD C1 integration of the solution and variational equation;
4. a strict interval Krawczyk inclusion.

The scripts report success only after the singular start, the complete flow,
the real G2 chamber, the full interval Jacobian and the strict Krawczyk
inclusion have all been verified.

## Requirements and reference configuration

- Julia 1.12;
- a C++17 compiler;
- CMake and Git for building CAPD;
- GMP and MPFR, including their development headers;
- `pkg-config`, required by the installed CAPD configuration script;
- CAPD compiled with multiprecision support and a working `capd-config`.

The Julia environment fixes `IntervalArithmetic` at version `1.0.10`.
A validation of both certificate scripts used the following configuration:

| Component | Version or setting |
| --- | --- |
| Julia | 1.12.6 |
| IntervalArithmetic | 1.0.10 |
| CAPD | 6.1.0, commit `37057e5899bccaa6cd5701057b7e6ea5e232ef8b` |
| CAPD backend | Native, with multiprecision enabled |
| C++ compiler | GCC 13.3.0 on Linux |
| GMP runtime | 6.3.0 |
| MPFR runtime | 4.2.1 |
| Julia and CAPD working precision | 256 bits |

### Install the pinned CAPD revision

For a new installation, use the specific revision below:

```bash
git clone https://github.com/CAPDGroup/CAPD.git CAPD
git -C CAPD checkout --detach 37057e5899bccaa6cd5701057b7e6ea5e232ef8b
cmake -S CAPD -B CAPD/build \
  -DCMAKE_INSTALL_PREFIX="$HOME/capd" \
  -DCAPD_ENABLE_MULTIPRECISION=ON
cmake --build CAPD/build -j
cmake --install CAPD/build

export CAPD_PREFIX="$HOME/capd"
export CAPD_CONFIG="$CAPD_PREFIX/bin/capd-config"
"$CAPD_CONFIG" --cflags
"$CAPD_CONFIG" --libs
```

The last two commands must succeed. If they report that `pkg-config` is
missing, install it before continuing. For an existing CAPD installation,
set `CAPD_CONFIG` to its configuration executable, or set `CAPD_PREFIX`
to its installation directory.

The Julia build script tries:

1. the executable named by `CAPD_CONFIG`;
2. `$CAPD_PREFIX/bin/capd-config`;
3. `capd-config` on `PATH`.

A working configuration script is required so that the compiler flags match
the installed CAPD backend. If `CAPD_PREFIX` is not set, it defaults to
`$HOME/capd`. On Apple silicon, GMP and MPFR are searched for under
`/opt/homebrew`; this can be changed by setting `GMP_PREFIX`.
The certificate flow rebuilds the C++ driver by default and prints the build
command, compiler version and executable SHA256 hash.

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
shooting residual. The floating-point right-hand-side evaluations check that
the square-root radicands are positive, and the final three chamber margins
are reported. The rigorous CAPD certificates check the three required chamber
signs throughout every validated integration step.
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
eight-by-eight matching Jacobian. At an exact match, its two time-derivative
columns agree: `dM/dT1 = Phi(q) = dM/dT2`. Thus `(dT1,dT2) = (1,-1)` is an
expected kernel direction. It moves the joining orbit along the same metric.
The unreduced matching derivative therefore cannot have rank eight at an
exact match. A computed rank of seven is consistent with this joining freedom;
it does not by itself establish a family of distinct metrics. Additional
kernel directions would require further analysis after fixing the joining
orbit and accounting for dependent matching conditions.

The matching residual must be read together with the reported start-extrapolation
error estimate. A residual smaller than that estimate does not establish the
same accuracy for the exact matching problem. Richardson extrapolation supplies
a numerical error estimate, rather than a rigorous interval bound. This program
finds candidates and tests the known squashed-S7 benchmark; it supplies no
rigorous existence, uniqueness or moduli certificate.

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

with radius `1e-10` in every coordinate. The default script applies
`radius_scale = 1/10000` to the base radii `1e-6` in `G2Equations.parameters`.
The resulting certificate establishes uniqueness in this smaller box.

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

Each certificate also prints the Julia and IntervalArithmetic versions,
working precision, hashes of `Project.toml`, `Manifest.toml` and the files in
`src/` and `deps/`, the CAPD build command, compiler version and executable hash.

### Save the final certificate logs

First commit the final source and scripts, including any newly added files.
Then run the commands below from the repository root. Use a shell supporting
`pipefail`, such as Bash or Zsh, so a failed Julia command remains a failed
command when its output is saved through `tee`:

```bash
mkdir -p logs
set -o pipefail
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

julia --startup-file=no --project=. scripts/test_capd_setup.jl 2>&1 | tee logs/capd_setup.log
julia --startup-file=no --project=. test/runtests.jl 2>&1 | tee logs/julia_tests.log
julia --startup-file=no --project=. scripts/certify_tau4_capd.jl 2>&1 | tee logs/tau4.log
julia --startup-file=no --project=. scripts/certify_tau23_u1_capd.jl 2>&1 | tee logs/tau23_u1.log

git rev-parse HEAD > logs/repository_commit.txt
git status --short > logs/working_tree_status.txt
git diff HEAD > logs/local_changes.patch
```

Keep the logs with that exact submission release. If any source or script
changes after a certificate run, rerun the affected certificate and update the
recorded revision. `working_tree_status.txt` should show no pending source or script changes;
untracked output logs are harmless. If an installed CAPD revision or the
GMP/MPFR versions differ from the reference configuration above, record the
versions used alongside the successful logs.

A certificate succeeds only if it exits successfully and prints the final
`SUCCESS: strict Krawczyk inclusion certified with CAPD C1 flow.` message.
A failed run is not an existence or uniqueness certificate.

## Main repository files

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
