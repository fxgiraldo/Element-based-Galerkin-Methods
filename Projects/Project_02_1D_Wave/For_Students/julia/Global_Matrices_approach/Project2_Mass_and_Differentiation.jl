using Plots

include("QuadraturePoints.jl")
include("element_matrices.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

#ipoints=1: Lobatto
#ipoints=2: Legendre
#ipoints=3: Chebyshev
#ipoints=4: Equi-spaced

#{{{ Main
function main()

    #Fix Interpolation and Integration Order
    N=1
    Q=N+1
    Np=N+1
    Nq=Q+1
    ipoints=1
    qpoints=1

    #Allocate Arrays
    M=zeros(DFloat,Np,Np)
    D=zeros(DFloat,Np,Np)

    #Select Interpolation Points
    if ipoints == 1
        (ξ,ω) = QuadraturePoints.lobatto_gauss(Np,DFloat) #Lobatto
    elseif ipoints == 2
        (ξ,ω) = QuadraturePoints.legendre_gauss(Np,DFloat) #Legendre
    elseif ipoints == 3
        (ξ,ω) = QuadraturePoints.chebyshev_gauss(Np,DFloat) #Chebyshev
    elseif ipoints == 4
        (ξ,ω) = QuadraturePoints.equispaced_points(Np,DFloat) #Equi-spaced
    end

    #Select Quadrature Points
    if qpoints == 1
        (ξq,ωq) = QuadraturePoints.lobatto_gauss(Nq,DFloat) #Lobatto
    elseif qpoints == 2
        (ξq,ωq) = QuadraturePoints.legendre_gauss(Nq,DFloat) #Legendre
    end

    #Construct Lagrange Polynomials
    (ψ,dψ) = QuadraturePoints.lagrange_basis(Np,Nq,ξ,ξq)

    (M,D) = element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)

    #Check correctness
    for j=1:Np, i=1:Np
        @show(i,j,M[i,j])
    end
    row_sum=zeros(DFloat,Np)
    for i=1:Np
        row_sum[i]=sum(D[i,:])
    end
    @show(row_sum)
    @show(sum(M))

    #Plot Interpolation
    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
