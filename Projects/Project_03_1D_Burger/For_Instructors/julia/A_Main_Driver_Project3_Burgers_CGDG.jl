#=
-------------------------------------------------------------------------------------------------------------
This file runs the 1D Burgers Equation using a unified CG/DG method using the AGGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.
Note that it only requires the storage of a global mass matrix, inside of TI_LSRK, 
the CREATE_RHS routine builds the RHS without using global matrices.

Written by F.X. Giraldo on July 6, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943

For time-integration, we use LSRK45 (4th order, 5-stage). 
For CG it runs OK until Time=0.5. For CG, it will blow up after but for DG it works for
quite a while longer.

The interpolation points used are the following:
#ipoints=1: Lobatto
#ipoints=2: Legendre
#ipoints=3: Chebyshev
#ipoints=4: Equi-spaced

The integration points used are the following:
#qpoints=1: Lobatto
#qpoints=2: Legendre
-------------------------------------------------------------------------------------------------------------
=#

using Plots

#Include all files used here
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

    #Fix Interpolation and Integration Order
    N=6
    Q=N
    Np=N+1
    Nq=Q+1
    ipoints=1
    qpoints=1
    space_method="CG" #time_final=0.4
    #space_method="DG" #time_final=0.4
    case=1 #1=exponential, 2=square wave
    Ne=21
    Npoin_dg=Ne*Np
    Npoin_cg=Ne*N + 1
    dt=DFloat(0.0001)
    time_final=DFloat(0.4)
    plot_movie=true
    Courant_max=0.5
    ntime = ceil(Int64, time_final / dt)
    @show(N,Q,Ne,Npoin_cg,Npoin_dg,dt,ntime,case)

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

    #Construct Local Mass and Differentiation Matrices
    (Me,De,Dwe,Fe) = element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)

    #Construct Grid and (i,e)->I map
    (coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg) = create_grid(Np,Npoin_cg,Npoin_dg,Ne,ξ,DFloat)

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

    #Construct Global Mass and Differentiation Matrices (only need Mass)
    (M,~,~,~) = global_matrices(Me,De,Dwe,Fe,intma,coord,Npoin,Ne,Np,periodicity,DFloat)
    
    #Compute Initial Condition
    time=DFloat(0.0)
    (qe) = exact_solution(coord,Npoin,time,case,DFloat)
    q0=copy(qe)
    dx_min=coord[2]-coord[1]

    #Recompute DT, NTIME
    umax=0
    for i=1:Npoin
        umax=max(umax, abs(q0[i]))
    end
    dt=Courant_max*dx_min/umax
    courant=umax*dt/dx_min
    ntime = ceil(Int64, time_final / dt)
    println("Courant = ",courant," dt = ",dt," ntime = ",ntime)

    #Time Integration
    (q0,time) = ti_LSRK!(q0,coord,M,Dwe,intma,periodicity,time,ntime,dt,space_method,plot_movie,DFloat)
    println(" ntime = ",ntime," time = ",time," qmax = ",maximum(q0)," qmin = ",minimum(q0))

    #Compute L2 Norm
    (qe) = exact_solution(coord,Npoin,time,case,DFloat)
    norms = compute_norms(q0,qe,Npoin,DFloat)
    println(" norm1 = ",norms[1]," norm2 = ",norms[2]," norm8 = ",norms[3])

    #Plot Interpolation
    qarray=zeros(DFloat,Npoin,2)
    qarray[:,1] .= q0
    qarray[:,2] .= qe
    plot_handle=plot(coord,qarray,xlabel="x",ylabel="q(x,t)",legend=true,lw=3,label=[space_method "Exact"],title=space_method)
    display(plot_handle)
#    savefig(plot_handle,"Project2.png")
    println("**Simulation Finished**") #output
end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
