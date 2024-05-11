#=
-------------------------------------------------------------------------------------------------------------
This file runs the 1D Wave Equation using a unified CG/DG method using the AGGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.
This code constructs global matrices, which although not the most efficient implementation, 
allows for the analysis of the resulting global matrices.

Written by F.X. Giraldo on July 5, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943

For time-integration, we use SSP methods of order 2 and 3; and LSRK45 (4th order, 5-stage). 

The interpolation points used are the following:
ipoints=1: Lobatto
ipoints=2: Legendre
ipoints=3: Chebyshev
ipoints=4: Equi-spaced

The integration points used are the following:
qpoints=1: Lobatto
qpoints=2: Legendre
-------------------------------------------------------------------------------------------------------------
=#

using Plots

include("QuadraturePoints.jl")
include("element_matrices.jl")
include("create_grid.jl")
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norms.jl")
include("ti_LSRK.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

#{{{ Main
function main()

    #------------Parameters you can change
    Ne=32
    N=4
    integration_type="exact" #exact or inexact
    ipoints=1 #1=LGL, 2=LG, 3=Cheby, 4=ESP
    qpoints=1 #1=LGL, 2=LG
    space_method="cg" #cg or dg
    time_final=1.0
    plot_movie=false
    #------------Parameters you can change
    
    #------------Parameters you CANNOT change
    Q=N
    if integration_type == "exact" && qpoints == 1 
        Q=N+1;
    end
    Np=N+1
    Nq=Q+1
    Npoin_dg=Ne*Np
    Npoin_cg=Ne*N + 1
    Courant_max=0.1
    case=1 #1=exponential, 2=square wave
    stages=5 #LSRK45 is 4th order with 5 stages
   #------------Parameters you CANNOT change

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

    #=--------------------Students Add These Files-----------------=#
    #=--------------------Students Add These Files-----------------=#
    #Construct Local Mass and Differentiation Matrices
    (Me,De,Dwe,Fe) = element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)
    #=--------------------Students Add These Files-----------------=#
    #=--------------------Students Add These Files-----------------=#

    #Construct Grid and (i,e)->I map
    (coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg) = create_grid(Np,Npoin_cg,Npoin_dg,Ne,ξ,DFloat)

    #CG or DG?
    if (space_method == "cg")
        Npoin=Npoin_cg
        coord=copy(coord_cg)
        intma=copy(intma_cg)
        periodicity=copy(periodicity_cg)
    elseif (space_method == "dg")
        Npoin=Npoin_dg
        coord=copy(coord_dg)
        intma=copy(intma_dg)
        periodicity=copy(periodicity_dg)
    end
    
    #=--------------------Students Add These Files-----------------=#
    #=--------------------Students Add These Files-----------------=#
    #Construct Global Mass and Differentiation Matrices
    (M,D,Dw,F) = global_matrices(Me,De,Dwe,Fe,intma,coord,Npoin,Ne,Np,periodicity,DFloat)
    #=--------------------Students Add These Files-----------------=#
    #=--------------------Students Add These Files-----------------=#
    Dhat=M\Dw
    Fhat=M\F

    #Compute Initial Condition
    time=DFloat(0.0)
    (qe,u) = exact_solution(coord,Npoin,time,case,DFloat)
    q0=copy(qe)
    dx_min=coord[2]-coord[1]

    #Recompute DT, NTIME
    dt=Courant_max*dx_min/u
    courant=u*dt/dx_min
    ntime = ceil(Int64, time_final / dt)
    #println("Courant = ",courant," dt = ",dt," ntime = ",ntime)
    @show(space_method,N,Q,Ne,Npoin,Npoin_cg,Npoin_dg,dt,ntime,time_final,courant,case,stages)

    #Time Integration
    (q0,time) = ti_LSRK!(q0,u,coord,Dhat,Fhat,intma,periodicity,time,time_final,ntime,dt,space_method,plot_movie,DFloat)
    println(" itime = ",ntime," dt = ",dt," time = ",time," qmax = ",maximum(q0)," qmin = ",minimum(q0))

    #Compute L2 Norm
    (qe,u) = exact_solution(coord,Npoin,time,case,DFloat)
    norms = compute_norms(q0,qe,Npoin,DFloat)
    println(" norm.L1 = ",norms[1]," norm.L2 = ",norms[2]," norm.L∞ = ",norms[3])

    #Plot Interpolation
    qarray=zeros(DFloat,Npoin,2)
    qarray[:,1] .= q0
    qarray[:,2] .= qe
    plot_handle=plot(coord,qarray,xlabel="x",ylabel="q(x,t)",legend=true,lw=3,label=[space_method "Exact"],title=space_method)
    display(plot_handle)
    #plot_solution(q0,coord,space_method,time)
#    savefig(plot_handle,"Project2.png")

    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
