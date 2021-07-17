#-------------------------------
#Allows reading in N
#-------------------------------
include("QuadraturePoints.jl")

#Some Constants
DFloat = Float16
N=parse(Int,ARGS[1])

#Some Output
@show (N)

# {{{ main
function main()

#Chebyshev Points
(ξ,ω) = QuadraturePoints.chebyshev_gauss(N+1)
println("Chebyshev Points ") #output
@show (ξ,ω)

#Legendre Points
(ξ,ω) = QuadraturePoints.legendre_gauss(N+1)
println("Legendre Points ") #output
@show (ξ,ω)

#Lobatto Points
(ξ,ω) = QuadraturePoints.lobatto_gauss(N+1)
println("Lobatto Points ") #output
@show (ξ,ω)

end

main()
