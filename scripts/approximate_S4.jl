#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "ApproximateShooting.jl"))
using .ApproximateShooting
using LinearAlgebra
using Printf

const DEFAULT_SEED = [0.04233, -0.06470, 0.00676, 0.76569]

function usage()
    println("Usage:")
    println("  julia --startup-file=no --project=. scripts/approximate_S4.jl")
    println("  julia --startup-file=no --project=. scripts/approximate_S4.jl a c nu T")
end

if any(arg -> arg in ("-h", "--help"), ARGS)
    usage()
    exit(0)
end

length(ARGS) in (0, 4) || (usage(); error("expected either zero or four arguments"))
seed = isempty(ARGS) ? copy(DEFAULT_SEED) : parse.(Float64, ARGS)

println("Non-rigorous numerical search for a zero of S4")
println("Initial guess (a,c,nu,T) = ", seed)
result = solve_approximate(:tau4, seed)

println("\nApproximate zero of S4:")
labels = ("a", "c", "nu", "T")
for (label, value) in zip(labels, result.parameters)
    @printf("  %-3s = %.17g\n", label, value)
end
println("S4 = ", result.residual)
@printf("infinity norm of S4 = %.6e\n", norm(result.residual, Inf))
println("singular start: delta = ", result.delta,
        ", Taylor order = ", result.start_order)
println("terminal chamber margins = ", result.chamber)
println("Newton iterations = ", result.iterations)
println("shooting-map evaluations = ", result.evaluations)
println("accepted steps in final flow = ", result.flow_steps)
println("\nThis computation is not a proof.")
println("Run scripts/certify_tau4_capd.jl for the rigorous certificate.")
