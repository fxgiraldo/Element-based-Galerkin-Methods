using Plots

include("QuadraturePoints.jl")
#include("element_matrices.jl")
include("create_grid.jl")
#include("global_matrices.jl")
include("exact_solution.jl")
#include("compute_norms.jl")
#include("time_integration.jl")

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
    N=3
    Q=N+1
    Np=N+1
    Nq=Q+1
    ipoints=1
    qpoints=1
    space_method="CG"
    #space_method="DG"
    #ti_type="SSP"
    ti_type="LSRK"
    case=2 #1=exponential, 2=square wave
    stages=3 #1,2, or 3 only
    if (ti_type == "LSRK")
        stages=5
    end
    Ne=40
    Npoin_dg=Ne*Np
    Npoin_cg=Ne*N + 1
    dt=DFloat(0.0001)
    time_final=DFloat(2.0)
    Courant_max=0.5
    ntime = ceil(Int64, time_final / dt)
    @show(N,Q,Ne,Npoin_cg,Npoin_dg,dt,ntime,case,ti_type,stages)

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
    (ψ,dψ) = QuadraturePoints.lagrange_basis(Np,Nq,ξ,ξq,DFloat)

    #Construct Grid and (i,e)->I map
    (coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg) = create_grid(Np,Npoin_cg,Npoin_dg,Ne,ξ,DFloat)

    #Construct Local Mass and Differentiation Matrices
#    (Me,De,Dwe,Fe) = element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)

    #CG or DG?
    if (space_method == "CG")
        Npoin=Npoin_cg
        coord=copy(coord_cg)
        intma=copy(intma_cg)
        periodicity=copy(periodicity_cg)
    elseif (space_method == "DG")
        Npoin=Npoin_dg
        coord=copy(coord_dg)
        intma=copy(intma_dg)
        periodicity=copy(periodicity_dg)
    end
    @show(space_method,Npoin,Ne,N,Q)

    #Construct Global Mass and Differentiation Matrices
#    (M,D,Dw,F) = global_matrices(Me,De,Dwe,Fe,intma,coord,Npoin,Ne,Np,periodicity,DFloat)

    #Compute Initial Condition
    time=DFloat(0.0)
    (qe,u) = exact_solution(coord,Npoin,time,case,DFloat)
    q0=copy(qe)
    dx_min=coord[2]-coord[1]

    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
