module ApproximateShooting

using LinearAlgebra
using IntervalArithmetic

# The numerical flow is ordinary Float64 arithmetic.  Correctly rounded
# intervals are used only to generate a reliable high-order initial value
# from the same regularized equations as the certificate.
IntervalArithmetic.configure(; rounding=:correct, matmul=:slow)

# Reuse the exact singular initial coefficients from the validated code.  The
# floating-point flow below is deliberately separate from the certificate.
include("ValidatedArithmetic.jl")
include("G2Equations.jl")

using .ValidatedArithmetic
using .G2Equations

export approximate_shooting, solve_approximate, chamber_margins
export vectorfield_float, high_order_start, integrate_float
export finite_difference_jacobian

const LAMBDA_FLOAT = 6.0 / sqrt(5.0)
const DEFAULT_ABSTOL = 1e-14
const DEFAULT_RELTOL = 1e-13

default_delta(kind::Symbol) =
    kind == :tau4 ? 1.0 / 8.0 :
    kind == :tau23 ? 1.0 / 1024.0 :
    kind == :tau23_u1 ? 1.0 / 1024.0 :
    error("Unknown shooting problem: $kind")

default_start_order(kind::Symbol) =
    kind == :tau4 ? 24 :
    kind == :tau23 ? 8 :
    kind == :tau23_u1 ? 8 :
    error("Unknown shooting problem: $kind")

function lift_float(u::AbstractVector, kind::Symbol)
    if kind in (:tau4, :tau23)
        length(u) == 8 || throw(DimensionMismatch("full state must have length 8"))
        return collect(u)
    elseif kind == :tau23_u1
        length(u) == 6 || throw(DimensionMismatch("tau23 U(1) state must have length 6"))
        q1, q2, q4, q5, q6, q8 = u
        return [q1, q2, q2, q4, q5, q6, q6, q8]
    end
    error("Unknown shooting problem: $kind")
end

select_float(q::AbstractVector, kind::Symbol) =
    kind in (:tau4, :tau23) ? collect(q) :
    kind == :tau23_u1 ? collect(q[[1, 2, 4, 5, 6, 8]]) :
    error("Unknown shooting problem: $kind")

function p_from_q_float(q::AbstractVector, lambda::Real)
    q1, q2, q3, q4, q5, q6, q7, q8 = q
    d1 = q2 + q7
    d2 = q3 + q6
    d3 = q4 + q5

    rad1 = -d1 * d2 / (lambda * d3)
    rad2 = -d1 * d3 / (lambda * d2)
    rad3 = -d2 * d3 / (lambda * d1)
    minimum((rad1, rad2, rad3)) > 0 ||
        throw(DomainError((rad1, rad2, rad3), "trajectory left the real G2 chamber"))

    return -sqrt(rad1), sqrt(rad2), sqrt(rad3)
end

"""Floating-point copy of the algebraically cancelled validated vector field."""
function vectorfield_float(u::AbstractVector, lambda::Real, kind::Symbol)
    q = lift_float(u, kind)
    q1, q2, q3, q4, q5, q6, q7, q8 = q
    p1, p2, p3 = p_from_q_float(q, lambda)

    a1 = q1 * q8
    a2 = q2 * q7
    a3 = q3 * q6
    a4 = q4 * q5
    sigma = a1 + a2 + a3 + a4
    ell = lambda / (2 * p1 * p2 * p3)

    f = [
        ell * (q1 * (2 * a1 - sigma) + 2 * q2 * q3 * q5),
        p3 - ell * (q2 * (2 * a2 - sigma) + 2 * q1 * q4 * q6),
        p2 - ell * (q3 * (2 * a3 - sigma) + 2 * q1 * q4 * q7),
        -p1 + ell * (q4 * (2 * a4 - sigma) + 2 * q2 * q3 * q8),
        p1 - ell * (q5 * (2 * a4 - sigma) + 2 * q1 * q6 * q7),
        -p2 + ell * (q6 * (2 * a3 - sigma) + 2 * q2 * q5 * q8),
        -p3 + ell * (q7 * (2 * a2 - sigma) + 2 * q3 * q5 * q8),
        -ell * (q8 * (2 * a1 - sigma) + 2 * q4 * q6 * q7),
    ]
    return select_float(f, kind)
end

function formal_state_coefficients(parameters::AbstractVector, kind::Symbol, order::Integer)
    order >= 2 || throw(ArgumentError("start order must be at least 2"))

    equation_kind = if kind in (:tau4, :tau23)
        start_parameters = parameters[1:3]
        :tau4
    elseif kind == :tau23_u1
        start_parameters = parameters[1:2]
        :tau23_u1
    else
        error("Unknown shooting problem: $kind")
    end

    p = [iv(x) for x in start_parameters]
    lambda = lambda_value()
    a, c = p[1:2]
    b, v = initial_pair(p, lambda, equation_kind)
    d = length(b)

    # For t*y'=G(t,y), the coefficient y_m satisfies
    # (m I-D_yG(0,y_0))y_m = [t^m]G(t,y_{<m}).
    linearization = jacobian(
        y -> regularized(y, a, c, iv(0), lambda, equation_kind), b,
    )
    coefficients = hcat(b, v, fill(iv(0), d, order - 1))

    for m in 2:order
        time = timejet(iv(0), m)
        y = [TJet(collect(coefficients[i, 1:m+1])) for i in 1:d]
        right_hand_side = regularized(y, a, c, time, lambda, equation_kind)
        forcing = [coefficient(right_hand_side[i], m) for i in 1:d]
        coefficients[:, m+1] = verified_solve(
            iv(m) * eye(d) - linearization, forcing,
        )
    end

    return coefficients
end

function high_order_start(parameters::AbstractVector, delta::Real,
                          kind::Symbol, order::Integer)
    delta > 0 || throw(ArgumentError("delta must be positive"))

    if kind == :tau4
        length(parameters) == 4 || throw(DimensionMismatch("tau4 parameters are (a,c,nu,T)"))
    elseif kind == :tau23
        length(parameters) == 4 || throw(DimensionMismatch("tau23 parameters are (a,c,nu,T)"))
    elseif kind == :tau23_u1
        length(parameters) == 3 || throw(DimensionMismatch("tau23 U(1) parameters are (a,c,T)"))
    else
        error("Unknown shooting problem: $kind")
    end

    a, c = parameters[1:2]
    a > 0 || throw(DomainError(a, "a must be positive"))
    c < 0 || throw(DomainError(c, "c must be negative"))

    coefficients = formal_state_coefficients(parameters, kind, order)
    equation_kind = kind == :tau23 ? :tau4 : kind
    delta_interval = iv(delta)
    y_delta_interval = [
        polynomial_value(TJet(collect(coefficients[i, :])), delta_interval)
        for i in axes(coefficients, 1)
    ]

    q_interval = q_from_y(
        y_delta_interval, iv(a), iv(c), delta_interval, equation_kind,
    )
    return Float64.(mid.(q_interval))
end

# One Dormand-Prince 5(4) step.  The fifth-order value advances the flow and
# the difference from the embedded fourth-order value controls the step size.
function dopri54_step(f, t, q, h)
    k1 = f(t, q)
    k2 = f(t + h / 5, q + h * (k1 / 5))
    k3 = f(t + 3h / 10, q + h * (3k1 / 40 + 9k2 / 40))
    k4 = f(t + 4h / 5, q + h * (44k1 / 45 - 56k2 / 15 + 32k3 / 9))
    k5 = f(t + 8h / 9,
           q + h * (19372k1 / 6561 - 25360k2 / 2187 +
                    64448k3 / 6561 - 212k4 / 729))
    k6 = f(t + h,
           q + h * (9017k1 / 3168 - 355k2 / 33 + 46732k3 / 5247 +
                    49k4 / 176 - 5103k5 / 18656))

    q5 = q + h * (35k1 / 384 + 500k3 / 1113 + 125k4 / 192 -
                  2187k5 / 6784 + 11k6 / 84)
    k7 = f(t + h, q5)
    q4 = q + h * (5179k1 / 57600 + 7571k3 / 16695 + 393k4 / 640 -
                  92097k5 / 339200 + 187k6 / 2100 + k7 / 40)

    return q5, q5 - q4
end

function integrate_float(q0::AbstractVector, t0::Real, t1::Real, kind::Symbol;
                         lambda=LAMBDA_FLOAT,
                         orientation=1.0,
                         abstol=DEFAULT_ABSTOL,
                         reltol=DEFAULT_RELTOL,
                         initial_step=1e-4, maxsteps=1_000_000)
    t1 > t0 || throw(ArgumentError("terminal time must exceed the start time"))
    q = Float64.(q0)
    t = Float64(t0)
    h = min(Float64(initial_step), Float64(t1 - t0))
    steps = 0
    orientation in (-1, 1, -1.0, 1.0) ||
        throw(ArgumentError("orientation must be +1 or -1"))
    rhs = (_, state) -> orientation .* vectorfield_float(state, lambda, kind)

    while t < t1
        steps += 1
        steps <= maxsteps || error("adaptive integrator exceeded maxsteps")
        h = min(h, t1 - t)

        trial, error_estimate = dopri54_step(rhs, t, q, h)
        all(isfinite, trial) || error("non-finite state produced by numerical flow")
        scale = abstol .+ reltol .* max.(abs.(q), abs.(trial))
        error_norm = maximum(abs.(error_estimate) ./ scale)

        if error_norm <= 1
            t += h
            q = trial
        end

        factor = error_norm == 0 ? 5.0 : clamp(0.9 * error_norm^(-0.2), 0.2, 5.0)
        h *= factor
        h > eps(t) || error("adaptive step size underflow")
    end

    return q, steps
end

function matching_float(q_reduced::AbstractVector, kind::Symbol)
    if kind == :tau4
        q = q_reduced
        return [q[1] - q[8], q[2] + q[4], q[3] - q[6], q[5] + q[7]]
    elseif kind == :tau23
        q = q_reduced
        return [q[1] - q[8], q[2] - q[7], q[3] - q[6], q[4] - q[5]]
    elseif kind == :tau23_u1
        q1, q2, q4, q5, q6, q8 = q_reduced
        return [q1 - q8, q2 - q6, q4 - q5]
    end
    error("Unknown shooting problem: $kind")
end

function chamber_margins(q_reduced::AbstractVector, kind::Symbol)
    q = lift_float(q_reduced, kind)
    return [q[2] + q[7], q[3] + q[6], -(q[4] + q[5])]
end

"""Evaluate the non-rigorous shooting map at one parameter vector."""
function approximate_shooting(kind::Symbol, parameters::AbstractVector;
                              delta=nothing, start_order=nothing,
                              abstol=DEFAULT_ABSTOL,
                              reltol=DEFAULT_RELTOL)
    delta_value = something(delta, default_delta(kind))
    order_value = something(start_order, default_start_order(kind))
    terminal_time = parameters[end]
    terminal_time > delta_value ||
        throw(DomainError(terminal_time, "terminal time must exceed delta"))

    q_delta = high_order_start(parameters, delta_value, kind, order_value)
    q_terminal, steps = integrate_float(
        q_delta, delta_value, terminal_time, kind;
        abstol=abstol, reltol=reltol,
    )
    residual = matching_float(q_terminal, kind)
    return (
        residual=residual,
        endpoint=q_terminal,
        chamber=chamber_margins(q_terminal, kind),
        steps=steps,
        delta=delta_value,
        start_order=order_value,
    )
end

function pack_parameters(parameters::AbstractVector, delta::Real, kind::Symbol)
    if kind in (:tau4, :tau23)
        a, c, nu, terminal_time = parameters
        return [log(a), log(-c), nu, log(terminal_time - delta)]
    elseif kind == :tau23_u1
        a, c, terminal_time = parameters
        return [log(a), log(-c), log(terminal_time - delta)]
    end
    error("Unknown shooting problem: $kind")
end

function unpack_parameters(x::AbstractVector, delta::Real, kind::Symbol)
    if kind in (:tau4, :tau23)
        return [exp(x[1]), -exp(x[2]), x[3], delta + exp(x[4])]
    elseif kind == :tau23_u1
        return [exp(x[1]), -exp(x[2]), delta + exp(x[3])]
    end
    error("Unknown shooting problem: $kind")
end

function finite_difference_jacobian(map, x; relative_step=2e-5)
    f0 = map(x)
    jacobian = Matrix{Float64}(undef, length(f0), length(x))
    for j in eachindex(x)
        h = relative_step * max(1.0, abs(x[j]))
        xp = copy(x)
        xm = copy(x)
        xp[j] += h
        xm[j] -= h
        jacobian[:, j] = (map(xp) - map(xm)) / (2h)
    end
    return jacobian
end

"""Locate an approximate zero by damped Newton iteration."""
function solve_approximate(kind::Symbol, seed::AbstractVector;
                           delta=nothing, start_order=nothing,
                           tolerance=1e-11, max_iterations=12,
                           relative_step=2e-5,
                           abstol=DEFAULT_ABSTOL,
                           reltol=DEFAULT_RELTOL,
                           verbose=true)
    delta_value = something(delta, default_delta(kind))
    order_value = something(start_order, default_start_order(kind))
    expected_length = kind in (:tau4, :tau23) ? 4 : kind == :tau23_u1 ? 3 :
        error("Unknown shooting problem: $kind")
    length(seed) == expected_length ||
        throw(DimensionMismatch("$kind requires $expected_length parameters"))
    x = pack_parameters(Float64.(seed), delta_value, kind)
    evaluations = 0

    function transformed_map(xvalue)
        evaluations += 1
        parameters = unpack_parameters(xvalue, delta_value, kind)
        return approximate_shooting(
            kind, parameters;
            delta=delta_value,
            start_order=order_value,
            abstol=abstol,
            reltol=reltol,
        ).residual
    end

    residual = transformed_map(x)
    for iteration in 0:max_iterations
        residual_norm = norm(residual, Inf)
        verbose && println("iteration ", iteration,
                           ": residual infinity norm = ", residual_norm)
        if residual_norm < tolerance
            parameters = unpack_parameters(x, delta_value, kind)
            shot = approximate_shooting(
                kind, parameters;
                delta=delta_value,
                start_order=order_value,
                abstol=abstol,
                reltol=reltol,
            )
            return (
                parameters=parameters,
                residual=shot.residual,
                endpoint=shot.endpoint,
                chamber=shot.chamber,
                flow_steps=shot.steps,
                evaluations=evaluations + 1,
                iterations=iteration,
                delta=delta_value,
                start_order=order_value,
            )
        end
        iteration == max_iterations && break

        jacobian = finite_difference_jacobian(
            transformed_map, x; relative_step=relative_step,
        )
        newton_step = -(jacobian \ residual)

        accepted = false
        damping = 1.0
        for _ in 1:14
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
                if !(exception isa DomainError || exception isa ArgumentError)
                    rethrow()
                end
            end
            damping /= 2
        end
        accepted || error("damped Newton iteration failed to decrease the residual")
    end

    error("Newton iteration did not converge; final residual infinity norm = $(norm(residual, Inf))")
end

end
