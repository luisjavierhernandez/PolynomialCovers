module PolynomialCovers


# =========================
# Basic structures
# =========================

include("Geometry.jl")

include("PolynomialUtils.jl")

include("PolynomialCriticalDisks.jl")

include("PreLifting.jl")

include("ApproximateLiftPath.jl")


include("FiberTransport.jl")



export PolynomialCover,
       coefficients,
       degree,
       value,

       Vertex,
       Path,
       Disk,
       distance,

       Fiber,
       critical_values_radius_origin,
       critical_values_disk_origin,
       fiber,
       PointState,
       PathState,

       basepoint,
       basepoint_initial_loop,

       safe_path_arc,
       avoid_points_path,

       transport!,
       transport_to!,
       cyclic_lift_solver
     

end # module
