#---------------------------------------------------------------------#
#This module stores everything needed to build different qudrature points
#in 1D.
#Written by F.X. Giraldo on 4/10/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
module QuadraturePoints

include("chebyshev_gauss.jl")
include("equispaced_points.jl")
include("legendre_poly.jl")
include("legendre_gauss.jl")
include("lobatto_gauss.jl")
#include("lagrange_basis.jl")

end
