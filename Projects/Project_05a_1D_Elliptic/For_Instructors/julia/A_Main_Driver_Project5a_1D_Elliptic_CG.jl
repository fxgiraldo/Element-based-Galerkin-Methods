#=
-------------------------------------------------------------------------------------------------------------
This code computes the 1D Poisson Equation using the CG method.

This is the solution to Project 5a in Element-based Galerkin Methods
The algorithms here are described in Algorithms 8.1 and 8.2, and DSSed
using Alg. 5.4. in F.X. Giraldo's Introduction to Element-based Galerkin 
Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.

Written by F.X. Giraldo on July 5, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943

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

#Include all files used here
include("QuadraturePoints.jl")
include("create_grid.jl")
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norms.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

function main()

    #-----------------------------Only Change these Input parameters---------------------------------#
    Ne=10 #Number of elements
    N=4 #Polynomial Order
    Q=N+1 #Quadrature Order
    ipoints=1 #interpolation point type
    qpoints=1 #integration point type
    space_method="CG"
    case=1 
    #-----------------------------Only Change these Input parameters---------------------------------#

    Np=N+1
    Nq=Q+1
    Npoin_dg=Ne*Np
    Npoin_cg=Ne*N + 1
    @show(N,Q,Ne,Npoin_cg,Npoin_dg,case)

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
    #-----------------------------Students Add their Function here--------------------------#
    #-----------------------------Students Add their Function here--------------------------#
    (M,L) = global_matrices(ψ,dψ,ωq,intma,coord,Npoin,Ne,Np,Nq,DFloat)
    #-----------------------------Students Add their Function here--------------------------#
    #-----------------------------Students Add their Function here--------------------------#
    
    #Compute Exact Solution
    (qe,qeₓ,fe) = exact_solution(coord,Npoin,case,DFloat)

    #Build RHS vector and apply Dirichlet BC
    R=M*fe
    R[1]=qe[1]
    R[Npoin]=qe[Npoin]

    #Compute CG Solution
    q0=-L\R

    #Compute L2 Norm
    norms = compute_norms(q0,qe,Npoin,DFloat)
    println(" norm1 = ",norms[1]," norm2 = ",norms[2]," norm∞ = ",norms[3])

    #Plot Interpolation
    p1 = plot(coord,qe,linecolor = :blue,label = "Analytic Solution",legend = :bottomright,legendfont=font(5))
    plot!(p1,coord,q0,linestyle = :dashdot,linecolor = :red,label = "$space_method Solution")
    title!(p1,"1D Poisson Problem for Ne= $Ne, N=$N, Q=$Q")
    xlabel!(p1,"x")
    ylabel!(p1,"q(x)")
    plot!(p1,xlims = (first(coord),last(coord)+.1),ylims = (-1,1))

    figure_1 = plot(p1)
    display(figure_1)

#    savefig(plot_handle,"Project2.png")
    println("**Simulation Finished**") #output
end

#----------------------------------#
# Run the main function
#----------------------------------#
main()
