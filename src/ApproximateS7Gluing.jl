module ApproximateS7Gluing

using LinearAlgebra
using ..ApproximateShooting

export approximate_s7_matching, solve_s7_matching, squashed_s7_seed
export xi_singular_start, pack_s7_parameters, unpack_s7_parameters

const VA = ApproximateShooting.ValidatedArithmetic
const LAMBDA_FLOAT = 6.0 / sqrt(5.0)
const DEFAULT_DELTA_ETA = 1.0 / 64.0
const DEFAULT_DELTA_XI = 3e-5
const DEFAULT_ETA_ORDER = 24

"""Known homogeneous squashed-S7 match, in the program's parameter order."""
function squashed_s7_seed()
    root5 = sqrt(5.0)
    return [
        root5 / 25, -3 * root5 / 5, root5 / 50, pi / 6,
        root5 / 25, -2 * root5 / 25, 19 * root5 / 25, pi / 6,
    ]
end

# ---------------------------------------------------------------------------
# Regular-singular start at the (1,-1) orbit
# ---------------------------------------------------------------------------

"""Closed-form leading coefficients y(0) at the (1,-1) orbit."""
function xi_initial_y0_guess(a::Real, b::Real, c::Real, lambda::Real)
    determinant = a * c - b^2
    difference = c - a
    determinant > 0 || throw(DomainError(determinant, "a*c-b^2 must be positive"))
    difference > 0 || throw(DomainError(difference, "c-a must be positive"))

    a_minus_c = a - c
    sqrt_determinant = sqrt(determinant)
    sqrt_difference = sqrt(difference)
    sqrt_product = sqrt(determinant * difference)
    sqrt_lambda = sqrt(lambda)
    lambda_squared = lambda^2
    lambda_three_halves = lambda * sqrt_lambda
    sqrt_difference_over_lambda = sqrt(difference / lambda)
    sqrt_lambda_difference_over_determinant =
        sqrt(lambda * difference / determinant)
    determinant_three_halves = determinant * sqrt_determinant
    difference_three_halves = difference * sqrt_difference
    denominator = 8 * determinant
    ac_minus_three_b_squared_over_two = a * c - 3 * b^2 / 2

    term1_r1 =
        -lambda_squared * (a^2 + a * c - 2 * b^2) * sqrt_determinant / 2 +
        (b^2 * sqrt_lambda_difference_over_determinant * lambda_three_halves +
         2 * a_minus_c * sqrt_lambda) * a * sqrt_difference
    inner_r1 = term1_r1 * sqrt_product -
        2 * a * lambda_squared * a_minus_c * sqrt_difference *
        ac_minus_three_b_squared_over_two
    r1 = -2 * inner_r1 / (difference_three_halves * denominator)

    r2 = sqrt_difference_over_lambda *
        (-sqrt_product + determinant * lambda_three_halves) / sqrt_product

    term_left_r3 =
        (-lambda_squared * sqrt_difference * a *
         ac_minus_three_b_squared_over_two / 2 +
         sqrt_lambda * determinant_three_halves) * sqrt_product
    term_bracket_r3 =
        -lambda_squared * (a^2 + a * c - 2 * b^2) * sqrt_determinant / 2 +
        (b^2 * sqrt_lambda_difference_over_determinant * lambda_three_halves -
         2 * a_minus_c * sqrt_lambda) * a * sqrt_difference
    inner_r3 = term_left_r3 - term_bracket_r3 * determinant / 4
    r3 = -8 * inner_r3 / (sqrt_product * sqrt_difference * denominator)

    r4 = -sqrt_difference_over_lambda *
        (determinant * lambda_three_halves + sqrt_product) / sqrt_product
    r5 = sqrt_difference_over_lambda *
        (-determinant * lambda_three_halves + sqrt_product) / sqrt_product

    term_left_r6 =
        (-c * lambda_squared * sqrt_difference *
         ac_minus_three_b_squared_over_two / 2 +
         sqrt_lambda * determinant_three_halves) * sqrt_product
    term_bracket_r6 =
        lambda_squared * (a * c - 2 * b^2 + c^2) * sqrt_determinant / 2 +
        (b^2 * sqrt_lambda_difference_over_determinant * lambda_three_halves +
         2 * a_minus_c * sqrt_lambda) * sqrt_difference * c
    inner_r6 = term_left_r6 - term_bracket_r6 * determinant / 4
    r6 = inner_r6 / (sqrt_product * sqrt_difference * determinant)

    r7 = sqrt_difference_over_lambda *
        (determinant * lambda_three_halves + sqrt_product) / sqrt_product

    term_bracket_r8 =
        lambda_squared * (a * c - 2 * b^2 + c^2) * sqrt_determinant / 2 +
        (b^2 * sqrt_lambda_difference_over_determinant * lambda_three_halves -
         2 * a_minus_c * sqrt_lambda) * sqrt_difference * c
    inner_r8 = 2 * term_bracket_r8 * sqrt_product -
        4 * lambda_squared * a_minus_c * sqrt_difference * c *
        ac_minus_three_b_squared_over_two
    r8 = inner_r8 / (difference_three_halves * denominator)

    return [r1, r2, r3, r4, r5, r6, r7, r8]
end

function xi_q_from_leading(y, a, b, c, time)
    time_squared = time^2
    return [
        a + time_squared * y[1],
        b + time * y[2],
        -a + time_squared * y[3],
        -b + time * y[4],
        b + time * y[5],
        c + time_squared * y[6],
        -b + time * y[7],
        -c + time_squared * y[8],
    ]
end

"""
Leading regular-singular start for xi.

The endpoint calculation applies Richardson extrapolation to starts at
`delta` and `delta/2`, cancelling the leading O(delta^2) start error.
"""
function xi_singular_start(parameters::AbstractVector, delta::Real;
                           lambda=LAMBDA_FLOAT)
    length(parameters) == 4 ||
        throw(DimensionMismatch("xi parameters are (a,b,c,T)"))
    a, b, c = parameters[1:3]
    a > 0 || throw(DomainError(a, "xi requires a>0"))
    c > a || throw(DomainError(c - a, "xi requires c>a"))
    a * c - b^2 > 0 ||
        throw(DomainError(a * c - b^2, "xi requires a*c-b^2>0"))
    delta > 0 || throw(DomainError(delta, "delta must be positive"))
    leading = xi_initial_y0_guess(a, b, c, lambda)
    return Float64.(xi_q_from_leading(leading, a, b, c, delta))
end

# ---------------------------------------------------------------------------
# Matching map and damped Newton solver
# ---------------------------------------------------------------------------

function validate_s7_parameters(parameters, delta_eta, delta_xi)
    length(parameters) == 8 || throw(DimensionMismatch(
        "parameters are (a_eta,c_eta,nu_eta,T1,a_xi,b_xi,c_xi,T2)",
    ))
    a_eta, c_eta, nu_eta, time_eta, a_xi, b_xi, c_xi, time_xi = parameters
    a_eta > 0 || throw(DomainError(a_eta, "eta requires a>0"))
    c_eta < 0 || throw(DomainError(c_eta, "eta requires c<0"))
    time_eta > delta_eta || throw(DomainError(time_eta, "T1 must exceed delta_eta"))
    a_xi > 0 || throw(DomainError(a_xi, "xi requires a>0"))
    c_xi > a_xi || throw(DomainError(c_xi - a_xi, "xi requires c>a"))
    a_xi * c_xi - b_xi^2 > 0 ||
        throw(DomainError(a_xi * c_xi - b_xi^2, "xi requires a*c-b^2>0"))
    time_xi > delta_xi || throw(DomainError(time_xi, "T2 must exceed delta_xi"))
    return nothing
end

"""
Evaluate eta(T1)-xi(T2) using independent local data at the two singular ends.

This is a floating-point candidate search, not an existence certificate.
"""
function approximate_s7_matching(parameters::AbstractVector;
                                 delta_eta=DEFAULT_DELTA_ETA,
                                 delta_xi=DEFAULT_DELTA_XI,
                                 eta_order=DEFAULT_ETA_ORDER,
                                 abstol=1e-14,
                                 reltol=1e-13)
    values = Float64.(parameters)
    validate_s7_parameters(values, delta_eta, delta_xi)
    eta_parameters = values[1:4]
    xi_parameters = values[5:8]

    eta_start = ApproximateShooting.high_order_start(
        eta_parameters, delta_eta, :tau4, eta_order,
    )
    eta_endpoint, eta_steps = ApproximateShooting.integrate_float(
        eta_start, delta_eta, eta_parameters[4], :tau4;
        orientation=1, abstol=abstol, reltol=reltol,
    )

    xi_coarse_start = xi_singular_start(
        xi_parameters, delta_xi; lambda=LAMBDA_FLOAT,
    )
    xi_coarse_endpoint, xi_coarse_steps = ApproximateShooting.integrate_float(
        xi_coarse_start, delta_xi, xi_parameters[4], :tau4;
        orientation=-1, abstol=abstol, reltol=reltol,
    )
    fine_delta = delta_xi / 2
    xi_fine_start = xi_singular_start(
        xi_parameters, fine_delta; lambda=LAMBDA_FLOAT,
    )
    xi_fine_endpoint, xi_fine_steps = ApproximateShooting.integrate_float(
        xi_fine_start, fine_delta, xi_parameters[4], :tau4;
        orientation=-1, abstol=abstol, reltol=reltol,
    )
    xi_endpoint = (4 * xi_fine_endpoint - xi_coarse_endpoint) / 3
    xi_extrapolation_error = norm(
        xi_fine_endpoint - xi_coarse_endpoint, Inf,
    ) / 3

    return (
        residual=eta_endpoint - xi_endpoint,
        eta_endpoint=eta_endpoint,
        xi_endpoint=xi_endpoint,
        eta_chamber=ApproximateShooting.chamber_margins(eta_endpoint, :tau4),
        xi_chamber=ApproximateShooting.chamber_margins(xi_endpoint, :tau4),
        eta_steps=eta_steps,
        xi_steps=xi_coarse_steps + xi_fine_steps,
        xi_extrapolation_error=xi_extrapolation_error,
        delta_eta=delta_eta,
        delta_xi=delta_xi,
        eta_order=eta_order,
    )
end

"""Map physical parameters to unconstrained Newton variables."""
function pack_s7_parameters(parameters::AbstractVector, delta_eta::Real,
                            delta_xi::Real)
    values = Float64.(parameters)
    validate_s7_parameters(values, delta_eta, delta_xi)
    a_eta, c_eta, nu_eta, time_eta, a_xi, b_xi, c_xi, time_xi = values
    ratio = b_xi / sqrt(a_xi * c_xi)
    abs(ratio) < 1 || throw(DomainError(ratio, "xi correlation must have magnitude <1"))
    return [
        log(a_eta), log(-c_eta), nu_eta, log(time_eta - delta_eta),
        log(a_xi), log(c_xi - a_xi), atanh(ratio), log(time_xi - delta_xi),
    ]
end

"""Map unconstrained Newton variables to admissible physical parameters."""
function unpack_s7_parameters(x::AbstractVector, delta_eta::Real,
                              delta_xi::Real)
    length(x) == 8 || throw(DimensionMismatch("S7 transformed state has length 8"))
    a_eta = exp(x[1])
    c_eta = -exp(x[2])
    nu_eta = x[3]
    time_eta = delta_eta + exp(x[4])
    a_xi = exp(x[5])
    c_xi = a_xi + exp(x[6])
    b_xi = sqrt(a_xi * c_xi) * tanh(x[7])
    time_xi = delta_xi + exp(x[8])
    return [a_eta, c_eta, nu_eta, time_eta, a_xi, b_xi, c_xi, time_xi]
end

function final_s7_result(x, residual, iterations, evaluations;
                         delta_eta, delta_xi, eta_order,
                         abstol, reltol, relative_step, rank_tolerance)
    parameters = unpack_s7_parameters(x, delta_eta, delta_xi)
    shot = approximate_s7_matching(
        parameters; delta_eta=delta_eta, delta_xi=delta_xi,
        eta_order=eta_order,
        abstol=abstol, reltol=reltol,
    )
    transformed_map = xvalue -> approximate_s7_matching(
        unpack_s7_parameters(xvalue, delta_eta, delta_xi);
        delta_eta=delta_eta, delta_xi=delta_xi,
        eta_order=eta_order,
        abstol=abstol, reltol=reltol,
    ).residual
    jacobian = ApproximateShooting.finite_difference_jacobian(
        transformed_map, x; relative_step=relative_step,
    )
    factorization = svd(jacobian)
    singular_values = factorization.S
    threshold = rank_tolerance * singular_values[1]
    numerical_rank = count(>(threshold), singular_values)
    condition_number = singular_values[1] / singular_values[end]
    kernel_basis = numerical_rank == length(singular_values) ?
        zeros(8, 0) : factorization.V[:, numerical_rank + 1:end]
    return (
        parameters=parameters,
        residual=shot.residual,
        eta_endpoint=shot.eta_endpoint,
        xi_endpoint=shot.xi_endpoint,
        eta_chamber=shot.eta_chamber,
        xi_chamber=shot.xi_chamber,
        eta_steps=shot.eta_steps,
        xi_steps=shot.xi_steps,
        xi_extrapolation_error=shot.xi_extrapolation_error,
        iterations=iterations,
        evaluations=evaluations + 18,
        delta_eta=delta_eta,
        delta_xi=delta_xi,
        eta_order=eta_order,
        jacobian=jacobian,
        singular_values=singular_values,
        numerical_rank=numerical_rank,
        condition_number=condition_number,
        kernel_basis=kernel_basis,
        rank_tolerance=rank_tolerance,
    )
end

"""Locate an approximate zero of the eight-component S7 matching map."""
function solve_s7_matching(seed::AbstractVector=squashed_s7_seed();
                           delta_eta=DEFAULT_DELTA_ETA,
                           delta_xi=DEFAULT_DELTA_XI,
                           eta_order=DEFAULT_ETA_ORDER,
                           tolerance=1e-9,
                           max_iterations=12,
                           relative_step=2e-5,
                           rank_tolerance=1e-7,
                           abstol=1e-14,
                           reltol=1e-13,
                           verbose=true)
    x = pack_s7_parameters(seed, delta_eta, delta_xi)
    evaluations = 0

    function transformed_map(xvalue)
        evaluations += 1
        return approximate_s7_matching(
            unpack_s7_parameters(xvalue, delta_eta, delta_xi);
            delta_eta=delta_eta, delta_xi=delta_xi,
            eta_order=eta_order,
            abstol=abstol, reltol=reltol,
        ).residual
    end

    residual = transformed_map(x)
    for iteration in 0:max_iterations
        residual_norm = norm(residual, Inf)
        verbose && println(
            "iteration ", iteration,
            ": gluing residual infinity norm = ", residual_norm,
        )
        if residual_norm < tolerance
            return final_s7_result(
                x, residual, iteration, evaluations;
                delta_eta=delta_eta, delta_xi=delta_xi,
                eta_order=eta_order,
                abstol=abstol, reltol=reltol,
                relative_step=relative_step, rank_tolerance=rank_tolerance,
            )
        end
        iteration == max_iterations && break

        jacobian = ApproximateShooting.finite_difference_jacobian(
            transformed_map, x; relative_step=relative_step,
        )
        factorization = svd(jacobian)
        cutoff = 1e-10 * factorization.S[1]
        newton_step = zeros(8)
        for index in eachindex(factorization.S)
            if factorization.S[index] > cutoff
                newton_step .-= (
                    dot(factorization.U[:, index], residual) /
                    factorization.S[index]
                ) .* factorization.V[:, index]
            end
        end
        accepted = false
        damping = 1.0

        for _ in 1:16
            trial_x = x + damping * newton_step
            try
                trial_residual = transformed_map(trial_x)
                if norm(trial_residual, Inf) < residual_norm
                    x = trial_x
                    residual = trial_residual
                    accepted = true
                    break
                end
            catch exception
                if !(exception isa DomainError || exception isa ArgumentError ||
                     exception isa VA.ValidationFailure)
                    rethrow()
                end
            end
            damping /= 2
        end
        accepted || error("damped Newton iteration failed to decrease the gluing residual")
    end

    error("S7 matching iteration did not converge; final infinity norm = $(norm(residual, Inf))")
end

end
