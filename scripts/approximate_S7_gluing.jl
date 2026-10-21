#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "ApproximateShooting.jl"))
include(joinpath(@__DIR__, "..", "src", "ApproximateS7Gluing.jl"))

using .ApproximateS7Gluing
using LinearAlgebra
using Printf

function usage()
    println("Usage:")
    println("  julia --startup-file=no --project=. scripts/approximate_S7_gluing.jl")
    println("  julia --startup-file=no --project=. scripts/approximate_S7_gluing.jl \\")
    println("      a_eta c_eta nu_eta T1 a_xi b_xi c_xi T2")
    println()
    println("With no arguments the program starts at the known squashed-S7 solution.")
    println("The two copies of a and c are independent shooting parameters.")
end

if any(argument -> argument in ("-h", "--help"), ARGS)
    usage()
    exit(0)
end

length(ARGS) in (0, 8) || (usage(); error("expected either zero or eight arguments"))
seed = isempty(ARGS) ? squashed_s7_seed() : parse.(Float64, ARGS)

println("Non-rigorous numerical search for eta(T1) = xi(T2)")
println("Parameter order:")
println("  (a_eta,c_eta,nu_eta,T1,a_xi,b_xi,c_xi,T2)")
println("Initial guess = ", seed)
println("The xi flow uses the reversed right-end equation q_s=-Phi(q).")

result = solve_s7_matching(seed)

println("\nApproximate S7 match:")
labels = ("a_eta", "c_eta", "nu_eta", "T1",
          "a_xi", "b_xi", "c_xi", "T2")
for (label, value) in zip(labels, result.parameters)
    @printf("  %-7s = %.17g\n", label, value)
end

println("\neta(T1)-xi(T2) = ", result.residual)
@printf("infinity norm of the matching residual = %.6e\n",
        norm(result.residual, Inf))
println("eta(T1) = ", result.eta_endpoint)
println("xi(T2)  = ", result.xi_endpoint)
println("eta terminal chamber margins = ", result.eta_chamber)
println("xi terminal chamber margins  = ", result.xi_chamber)
println("eta start: delta = ", result.delta_eta,
        ", Taylor order = ", result.eta_order)
println("xi starts: delta = ", result.delta_xi,
        " and ", result.delta_xi / 2,
        " (leading regular-singular expansion with Richardson extrapolation)")
@printf("xi start-extrapolation error estimate = %.6e\n",
        result.xi_extrapolation_error)
println("accepted eta/xi flow steps = ",
        (result.eta_steps, result.xi_steps))
println("Newton iterations = ", result.iterations)
println("matching-map evaluations = ", result.evaluations)

println("\nSingular values of the transformed 8-by-8 matching Jacobian:")
println(result.singular_values)
@printf("estimated condition number = %.6e\n", result.condition_number)
println("numerical rank = ", result.numerical_rank,
        " (relative threshold ", result.rank_tolerance, ")")
println("At an exact match, moving the joining orbit gives an expected kernel.")
println("In the time coordinates this direction is (dT1,dT2)=(1,-1).")
if result.numerical_rank == 8
    println("An exact match cannot have rank eight in this unreduced formulation.")
    println("The reported numerical rank does not establish local isolation.")
else
    println("Numerical kernel dimension = ", 8 - result.numerical_rank)
    println("Kernel basis in the transformed parameter coordinates:")
    display(result.kernel_basis)
    println("The joining freedom alone accounts for one kernel direction.")
end
println("These numerical diagnostics do not determine geometric rigidity or moduli.")

if norm(result.residual, Inf) <= result.xi_extrapolation_error
    println("\nThe residual is no larger than the estimated xi start error.")
    println("It therefore does not establish that accuracy for the exact matching problem.")
end
println("The Richardson start error is an estimate, not a rigorous interval bound.")

println("\nThis floating-point computation locates candidates only; it is not a proof.")
