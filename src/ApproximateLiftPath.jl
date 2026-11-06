############################################################
# ApproximateLiftPath.jl
#
# Path lifting algorithms for polynomial covers.
#
# These routines implement a predictor–corrector strategy
# for lifting paths in the codomain while avoiding the
# critical value disk.
############################################################

"""
    aux_lift_predictor_trace!(P, state, ypath; S, max_subdiv=10)

Predictor phase for path lifting with trajectory recording.

This routine performs a sequence of predictor steps to lift
the polygonal path `ypath` under the polynomial cover `P`.
The method follows a linearized lifting rule
    x_{i+1} = x_i + (y_{i+1} - w_i)/P'(x_i),
while ensuring that the lifted segment remains outside
the critical value disk of radius `3S/2`.

If the safety condition fails, the step is subdivided
until a safe predictor step is found.

# Arguments
- `P::PolynomialCover` : polynomial cover
- `state::PathState` : current lifting state
- `ypath::Path` : polygonal path in the codomain
- `S::Real` : critical value exclusion radius

# Keyword Arguments
- `max_subdiv` : maximum allowed subdivisions of a segment

# Returns
`(state, success)`

- `state` : updated lifting state
- `success::Bool` : `true` if the predictor phase succeeded,
  `false` if subdivision exceeded `max_subdiv`.
"""
function aux_lift_predictor_trace!(
    P::PolynomialCover,
    state::PathState{T},
    ypath::Path{T};
    S::Real,
    max_subdiv::Int = 10,
) where {T<:Real}

    x = state.x
    w = state.w

    verts = ypath.vertices
    n = length(verts)

    i = 1
    subdiv = 0

    while i < n

        y_i   = verts[i].z
        y_ip1 = verts[i+1].z

        dP = derivative(P, x)
        abs(dP) < eps(T) && error("Near critical point")

        x_new = x + (y_ip1 - w)/dP
        w_new = value(P, x_new)

        if segment_intersects_disk(w, w_new, 3S/2) ||
        segment_intersects_disk(w_new, y_ip1, 3S/2)

            subdiv += 1
            subdiv > max_subdiv && return state, false

            y_mid = (y_i + y_ip1)/2

            x = x + (y_mid - w)/dP
            w = value(P, x)

            push!(state.pathx, x)
            push!(state.pathy, w)

            continue
        end

        x = x_new
        w = w_new

        push!(state.pathx, x)
        push!(state.pathy, w)

        i += 1
    end

    state.x = x
    state.w = w

    return state, true
end

"""
    refine_endpoint_newton_trace!(P, state, w_target; tol=1e-9, max_iter=100)

Newton refinement of the endpoint of a lifted path.

After the predictor phase, the endpoint of the lifted path
is refined using Newton's method to solve
    P(x) = w_target.

within the prescribed tolerance.

The intermediate Newton iterates are recorded in `state.pathx`
while the codomain value remains fixed at `w_target`.

# Arguments
- `P::PolynomialCover`
- `state::PathState`
- `w_target::Complex`

# Keyword Arguments
- `tol` : relative Newton tolerance
- `max_iter` : maximum Newton iterations

# Returns
`x` : refined endpoint in the domain.
"""
function refine_endpoint_newton_trace!(
    P::PolynomialCover,
    state::PathState{T},
    w_target::Complex{T};
    tol::Real = 1e-9,
    max_iter::Int = 100,
) where {T<:Real}

    x = state.x

    for _ in 1:max_iter

        Px = value(P, x)
        err = Px - w_target

        if abs(err) ≤ tol*(1 + abs(w_target))
            state.x = x
            state.w = w_target

            push!(state.pathx, x)
            push!(state.pathy, w_target)

            return x
        end

        dPx = derivative(P, x)
       abs(dPx) ≤ sqrt(eps(T)) * max(1, abs(x)) &&
            error("Derivative nearly vanished during Newton refinement")

        x -= err/dPx

        push!(state.pathx, x)
        push!(state.pathy, w_target)

    end

    error("Newton refinement did not converge in $max_iter iterations")
end

"""
    lift_path_trace!(P, state, ypath; S, newton_tol=1e-9)

Lift a polygonal path under a polynomial cover while
recording the full lifting trajectory.

This routine combines a predictor phase with a Newton
endpoint correction in order to compute a lifted path
starting from the current state `(x,w)`.

# Arguments
- `P::PolynomialCover`
- `state::PathState`
- `ypath::Path` : polygonal path in the codomain
- `S::Real` : critical value exclusion radius

# Keyword Arguments
- `newton_tol` : tolerance for Newton endpoint refinement

# Returns
`(state, success)`

- `state` : updated lifting state
- `success::Bool` : indicates whether lifting succeeded

# Notes
This routine stores the entire lifting trajectory and is
therefore well suited for **visualization and educational
examples**.

For faster lifting without storing intermediate points
see `lift_path_fast!`.
# Algorithm
The lifting is computed using a predictor–corrector strategy:

1. A homotopy-safe predictor step lifts each segment of the path.
2. A Newton correction ensures that the final point lies in the
   correct fiber of the polynomial cover.
"""
function lift_path_trace!(
    P::PolynomialCover,
    state::PathState{T},
    ypath::Path{T};
    S::Real,
    newton_tol::Real = 1e-9,
) where {T<:Real}

    state_pred, ok = aux_lift_predictor_trace!(P, state, ypath; S=S)

    ok || return state_pred, false

    w_target = ypath.vertices[end].z

    refine_endpoint_newton_trace!(
        P,
        state_pred,
        w_target;
        tol=newton_tol,
    )

    return state_pred, true
end

function aux_lift_predictor!(
    P::PolynomialCover,
    state::PointState{T},
    ypath::Path{T};
    S::Real,
    max_subdiv::Int = 10,
) where {T<:Real}

    x = state.x
    w = state.w

    verts = ypath.vertices
    n = length(verts)

    i = 1
    subdiv = 0

    while i < n

        y_i   = verts[i].z
        y_ip1 = verts[i+1].z

        # predictor step
        dP = derivative(P, x)
        abs(dP) < eps(T) && error("Near critical point")

        x_new = x + (y_ip1 - w)/dP
        w_new = value(P, x_new)

        # homotopy safety check
        if segment_intersects_disk(w_new, y_ip1, 3S/2)

            subdiv += 1
            subdiv > max_subdiv && return state, false

            # halve step in codomain
            y_mid = (y_i + y_ip1)/2

            # do predictor only for half-step
            x = x + (y_mid - w)/dP
            w = value(P, x)

            continue
        end

        # accept step
        x = x_new
        w = w_new
        i += 1
    end

    state.x = x
    state.w = w

    return state, true
end

function refine_endpoint_newton!(
    P::PolynomialCover,
    state::PointState{T},
    w_target::Complex{T};
    tol::Real = 1e-9,
    max_iter::Int = 100,
) where {T<:Real}

    x = state.x

    for _ in 1:max_iter
        Px = value(P, x)
        err = Px - w_target

        if abs(err) ≤ tol*(1 + abs(w_target))
            state.x = x
            state.w = w_target
            return x
        end

        dPx = derivative(P, x)
        abs(dPx) < sqrt(eps(T)) &&
            error("Derivative nearly vanished during Newton refinement")

        x -= err/dPx
    end

    error("Newton refinement did not converge in $max_iter iterations")
end

function lift_path_fast!(
    P::PolynomialCover,
    state::PointState{T},
    ypath::Path{T};
    S::Real,
    newton_tol::Real = 1e-9,
) where {T<:Real}

    # --- predictor lifting (homotopy-safe) ---
    state_pred, ok = aux_lift_predictor!(P, state, ypath; S=S)

    ok || return state_pred, false

    # --- endpoint Newton refinement ---
    w_target = ypath.vertices[end].z
    refine_endpoint_newton!(
        P,
        state_pred,
        w_target;
        tol=newton_tol,
    )

    return state_pred, true
end

