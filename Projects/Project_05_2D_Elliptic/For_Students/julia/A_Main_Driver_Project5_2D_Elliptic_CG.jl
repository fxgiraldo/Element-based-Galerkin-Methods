#=
-------------------------------------------------------------------------------------------------------------
This file contains the Student template for solving Project 5: the 2D Poisson Equation using the CG method with GGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.

The approached follows Algorithm 12.18 in the book.  

Written by F.X. Giraldo on July 10, 2021. Updated on May 30, 2024.
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
using LinearAlgebra

include("QuadraturePoints.jl")
include("create_grid.jl")
include("compute_metrics.jl")
include("element_matrices.jl")
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norm.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

#Main Function
function main()

    #-----------------------------Only Change these Input parameters---------------------------------#
    Nel=4
    N=16
    integration_type="exact" #exact or inexact
    c=1.0 #Constant in Exact solution
    plot_grid=true
    plot_solution=true
    warp_grid=true
    #-----------------------------Only Change these Input parameters---------------------------------#
    
    ipoints=1 #keep it at 1 (LGL)
    qpoints=1 #keep it at 1 (LGL) but can also do 2 (LG)
    #space_method="CG" #only capable of doing CG
    Q=N
    if integration_type == "exact"
        Q=N+1;
    end
    Nex=Nel
    Ney=Nel
    Np=N+1
    Nq=Q+1
    Npts=Np^2
    Nqs=Nq^2
    Ne=Nex*Ney
    Nx=Nex*N+1
    Ny=Ney*N+1
    Npoin=Nx*Ny
    Nboun=2*Nx + 2*(Ny-2)
    @show(N,Q,Nel,Ne,Npoin,Nboun,c)

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
    (coord,intma,iboun) = create_grid(Np,Npoin,Ne,Nboun,Nex,Ney,ξ,plot_grid,warp_grid,DFloat)

    #Construct Metric Terms
    (ξ_x,ξ_y,η_x,η_y,jac) = compute_metrics(coord,intma,ψ,dψ,Ne,Np,Nq,DFloat)

    @time begin
    #-----------------------------Students Add their Functions Here--------------------------#
    #-----------------------------Students Add their Functions Here--------------------------#
    #Can do Element Matrices followed by Global:
    #Construct Element Matrices
    #(Me,Le) = element_matrices(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,Ne,Np,Nq,DFloat)
    #Construct Global Matrices
    #(M,L) = global_matrices(Me,Le,intma,Ne,Np,Npoin,DFloat)
    #
    #OR, do the Global directly - Much Faster Code!!
    #
    #Construct Element and Global Matrices all at once
    (M,L) = global_matrices(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)
    #-----------------------------Students Add their Functions Here--------------------------#
    #-----------------------------Students Add their Functions Here--------------------------#

    #Compute Exact Solution
    (qe,fe) = exact_solution(coord,Npoin,c,DFloat)

    #Apply Dirichlet BCs: Impose exact solution
    R=M*fe
    for i=1:Nboun
        I=iboun[i]
        L[I,:] = zeros(DFloat,Npoin)
        L[I,I] = DFloat(1)
        R[I]=qe[I]
    end#i

    #Solve Matrix problem
    q=L\R
    end

    #Compute L2 Norm
    norm2 = compute_norm(q,qe,Npoin,DFloat)
    println("norm.L2 = ",norm2)

    #Plot Solution
    if (plot_solution)
        x=reshape(coord[1,:],(Nx,Ny))
        y=reshape(coord[2,:],(Nx,Ny))
        qq=reshape(q,(Nx,Ny))
        xx=x[:,1]
        yy=y[1,:]

        p1 = contourf(xx,yy,qq)
        title!(p1,"Numerical Solution")
        xlabel!(p1,"x")
        ylabel!(p1,"y")
 
        qq=reshape(qe,(Nx,Ny))
        p2 = contourf(xx,yy,qq)
        title!(p2,"Exact Solution")
        xlabel!(p2,"x")
        ylabel!(p2,"y")

        figure_1 = plot(p1,p2,layout = (1,2))
        display(figure_1)
    end

    println("*Done*") #output

end

#----------------------------------#
# Run the main function
#----------------------------------#
main()
