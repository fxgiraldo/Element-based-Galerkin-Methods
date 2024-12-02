#=
-------------------------------------------------------------------------------------------------------------
This file runs the 1D Burgers Equation using a unified CG/DG method using the AGGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.
Note that it only requires the storage of a global mass matrix, inside of TI_LSRK, 
the CREATE_RHS routine builds the RHS without using global matrices.

The strategy followed is very similar to that described in Alg. 7.6 except that we do not build the 
Flux matrix explicitly. This is all done in CREATE_RHS.jl.

Written by F.X. Giraldo on July 6, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943

For time-integration, we use LSRK45 (4th order, 5-stage). 
For CG it runs OK until Time=0.5. For CG, it will blow up after but for DG it works for
quite a while longer.

Lobatto (LGL) interpolation and quadrature points are used. 
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

function main()

    #-----------------------------Only Change these Input parameters---------------------------------#
    Ne=21
    N=6
    Q=N
    space_method="DG" #time_final=0.4 - CG or DG
    cgdg_method="strong" #strong or weak. Entropy-stability only works with Strong
    flux_method="ES" #ES=entropy stable, NES=not ES
    case=1 #1=exponential, 2=square wave
    time_final=DFloat(1.0)
    plot_movie=true
    Courant_max=DFloat(0.5)
    #-----------------------------Only Change these Input parameters---------------------------------#

    Np=N+1
    Nq=Q+1
    Npoin_dg=Ne*Np
    Npoin_cg=Ne*N + 1
    dt=DFloat(0.0001)
    ntime = ceil(Int64, time_final / dt)
    @show(N,Q,Ne,Npoin_cg,Npoin_dg,dt,ntime,case)

    #Select Interpolation Points
    (ξ,ω) = QuadraturePoints.lobatto_gauss(Np,DFloat) #Lobatto
    
    #Select Quadrature Points
    (ξq,ωq) = QuadraturePoints.lobatto_gauss(Nq,DFloat) #Lobatto
    
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
    @show(space_method,cgdg_method,flux_method,Npoin)

    #Construct Global Mass and Differentiation Matrices (only need Mass)
    (M,~,~,~) = global_matrices(Me,De,Dwe,Fe,intma,coord,Npoin,Ne,Np,periodicity,DFloat)
    
    #Compute Initial Condition
    time=DFloat(0.0)
    (qe) = exact_solution(coord,Npoin,time,case,DFloat)
    q0=copy(qe)
    dx_min=coord[2]-coord[1]

    #Recompute DT, NTIME
    (dt,courant)=compute_courant(q0,Courant_max,dx_min,DFloat) 
    ntime = ceil(Int64, time_final / dt)
    println("Courant = ",courant," dt = ",dt," ntime = ",ntime)

    #Time Integration
    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#
    (q0,time,itime) = ti_LSRK!(q0,coord,M,De,Dwe,intma,periodicity,time,time_final,Courant_max,space_method,cgdg_method,flux_method,plot_movie,DFloat)
    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#
    
    #Print Time and Extrema
    println(" itime = ",itime," time = ",time," Courant = ",courant," dt = ",dt," qmax = ",maximum(abs,q0))

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

#----------------------------------#
# Run the main function
#----------------------------------#
main()
