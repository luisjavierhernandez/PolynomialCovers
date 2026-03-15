

@testset "Path lifting with trace" begin

    # Polinomio P(x) = x^2
    # coeficientes en grado decreciente
    P = PolynomialCover([1.0,0.0,0.0])

    # punto inicial
    x0 = 2.0 + 0im
    w0 = value(P, x0)

    state = PathState(x0, w0)

    # camino en el codominio
    yverts = [
        Vertex(w0),
        Vertex(w0 + 0.5),
        Vertex(w0 + 1.0),
    ]

    ypath = Path(yverts)

    # radio crítico
    S = 0.5

    # ---------- predictor ----------
    state_pred, ok = PolynomialCovers.aux_lift_predictor_trace!(P, state, ypath; S=S)

    @test ok == true
    @test length(state_pred.pathx) ≥ 2
    @test length(state_pred.pathx) == length(state_pred.pathy)

    # el último punto debe aproximar la imagen objetivo
    @test abs(value(P, state_pred.x) - state_pred.w) < 1e-8


    # ---------- Newton refinement ----------

    w_target = ypath.vertices[end].z

    x_ref = PolynomialCovers.refine_endpoint_newton_trace!(
        P,
        state_pred,
        w_target;
        tol=1e-10
    )

    @test abs(value(P, x_ref) - w_target) < 1e-10

    # se registraron iteraciones
    @test length(state_pred.pathx) ≥ 3


    # ---------- lifting completo ----------

    state2 = PathState(x0, w0)

    state_final, ok2 = PolynomialCovers.lift_path_trace!(
        P,
        state2,
        ypath;
        S=S,
        newton_tol=1e-10
    )

    @test ok2 == true

    w_target = ypath.vertices[end].z

    # endpoint correcto
    @test abs(value(P, state_final.x) - w_target) < 1e-10

    # trayectoria almacenada
    @test length(state_final.pathx) == length(state_final.pathy)
    @test length(state_final.pathx) ≥ 3

end