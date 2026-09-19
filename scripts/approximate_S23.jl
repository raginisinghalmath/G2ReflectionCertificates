#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "ApproximateShooting.jl"))
using .ApproximateShooting
using LinearAlgebra
using Printf

# The deliberately nonzero value of nu lets the full shooting problem reveal
# the enhanced U(1) symmetry instead of imposing it in advance.
const DEFAULT_SEED = [0.02588, -0.07764, 0.005, 0.85780]

function usage()
    println("Usage:")
    println("  julia --startup-file=no --project=. scripts/approximate_S23.jl")
    println("  julia --startup-file=no --project=. scripts/approximate_S23.jl a c nu T")
    println("The full four-variable shooting problem is solved; nu is not fixed.")
end

if any(arg -> arg in ("-h", "--help"), ARGS)
    usage()
    exit(0)
end

length(ARGS) in (0, 4) || (usage(); error("expected either zero or four arguments"))
seed = isempty(ARGS) ? copy(DEFAULT_SEED) : parse.(Float64, ARGS)

println("Non-rigorous numerical search for a zero of the full S23 map")
println("Initial guess (a,c,nu,T) = ", seed)
result = solve_approximate(:tau23, seed)

println("\nApproximate zero of S23:")
labels = ("a", "c", "nu", "T")
for (label, value) in zip(labels, result.parameters)
    @printf("  %-3s = %.17g\n", label, value)
end
println("S23 = ", result.residual)
@printf("infinity norm of S23 = %.6e\n", norm(result.residual, Inf))
println("singular start: delta = ", result.delta,
        ", Taylor order = ", result.start_order)
println("terminal chamber margins = ", result.chamber)
println("Newton iterations = ", result.iterations)
println("shooting-map evaluations = ", result.evaluations)
println("accepted steps in final flow = ", result.flow_steps)
@printf("\nThe computation finds nu approximately %.6e.\n", result.parameters[3])
println("This numerical observation motivates the exact U(1) reduction used in the certificate.")
println("\nThis computation is not a proof.")
println("Run scripts/certify_tau23_u1_capd.jl for the rigorous certificate.")
