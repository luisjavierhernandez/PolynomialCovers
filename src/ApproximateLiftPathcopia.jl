########################################
# ApproximateLiftPath.jl
########################################





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

###Versiones anteriores de lift_path! (con traza explícita)


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

            # predictor half-step
            x = x + (y_mid - w)/dP
            w = value(P, x)

            # record intermediate point
            push!(state.pathx, x)
            push!(state.pathy, w)

            continue
        end

        # accept step
        x = x_new
        w = w_new

        # record step
        push!(state.pathx, x)
        push!(state.pathy, w)

        i += 1
    end

    state.x = x
    state.w = w

    return state, true
end

function refine_endpoint_newton_trace!(
    P::PolynomialCover,
    state::PathState{T},
    w_target::Complex;
    tol::Real = 1e-9,
    max_iter::Int = 100,
) where {T<:Real}

    w_target = Complex{T}(w_target)
#function refine_endpoint_newton_trace!(
 #   P::PolynomialCover,
  #  state::PathState{T},
  #  w_target::Complex{T};
  #  tol::Real = 1e-9,
  #  max_iter::Int = 100,
#) where {T<:Real}
    x = state.x

    for _ in 1:max_iter

        Px = value(P, x)
        err = Px - w_target

        if abs(err) ≤ tol*(1 + abs(w_target))
            state.x = x
            state.w = w_target

            # record final point
            push!(state.pathx, x)
            push!(state.pathy, w_target)

            return x
        end

        dPx = derivative(P, x)
        abs(dPx) < sqrt(eps(T)) &&
            error("Derivative nearly vanished during Newton refinement")

        x -= err/dPx

        # record Newton step
        push!(state.pathx, x)
        push!(state.pathy, w_target)

    end

    error("Newton refinement did not converge in $max_iter iterations")
end

function lift_path_trace!(
    P::PolynomialCover,
    state::PathState{T},
    ypath::Path{T};
    S::Real,
    newton_tol::Real = 1e-9,
) where {T<:Real}

    # --- predictor lifting (homotopy-safe) ---
    state_pred, ok = aux_lift_predictor_trace!(P, state, ypath; S=S)

    ok || return state_pred, false

    # --- endpoint Newton refinement ---
    w_target = ypath.vertices[end].z

    refine_endpoint_newton_trace!(
        P,
        state_pred,
        w_target;
        tol = newton_tol,
    )

    return state_pred, true
end