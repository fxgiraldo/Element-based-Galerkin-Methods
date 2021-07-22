#=
-------------------------------------------------------------------------------------------------------------
This file contains the solution to Project 7: the 3D Poisson Equation using the CG method with GGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.

The approached follows Algorithm 12.18 in the book but extended to 3D.

Written by F.X. Giraldo on July 21, 2021.
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
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norm.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

#Main Function
function main()

    #-----------------------------Only Change these Input parameters---------------------------------#
    Nex=2
    Ney=2
    Nez=2
    N=4
    Q=N+1
    ipoints=1
    qpoints=1
    space_method="CG"
    c=DFloat(1) #Constant in Exact solution
    plot_grid=true
    plot_solution=false
    warp_grid=true
    #-----------------------------Only Change these Input parameters---------------------------------#

    Np=N+1
    Nq=Q+1
    Ne=Nex*Ney*Nez
    Nx=Nex*N+1
    Ny=Ney*N+1
    Nz=Nez*N+1
    Npoin=Nx*Ny*Nz
    Nboun=2*Nx*Nz + 2*(Ny-2)*Nz + 2*(Nx-2)*(Ny-2)
    @show(N,Q,Ne,Npoin,Nboun,c)

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
    (coord,intma,iboun) = create_grid(Np,Npoin,Ne,Nboun,Nex,Ney,Nez,ξ,plot_grid,warp_grid,DFloat)

    #Construct Metric Terms
    (ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac) = compute_metrics(coord,intma,ψ,dψ,Ne,Np,Nq,DFloat)

    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#
    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#
    #Construct Global Matrices
    (M,L) = global_matrices(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)
    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#
    #-----------------------------Students Add their Functions inside ti_LSRK--------------------------#

    #Compute Exact Solution
    (qe,fe) = exact_solution(coord,Npoin,c,DFloat)

    #Apply Dirichlet BCs: Impose exact solution
    R=M*fe
    for i=1:Nboun
        I=iboun[i]
        L[I,:] = zeros(DFloat,Npoin)
        L[I,I] = DFloat(1.0)
        R[I]=qe[I]
    end#i

    #Solve Matrix problem
    q=L\R

    #Compute L2 Norm
    norm2 = compute_norm(q,qe,Npoin,DFloat)
    println("norm.L2 = ",norm2)

    #Plot Solution
    if (plot_solution)
        x=reshape(coord[1,:],(Nx,Ny,Nz))
        y=reshape(coord[2,:],(Nx,Ny,Nz))
        qq=reshape(q,(Nx,Ny,Nz))
        xx=x[:,1,1]
        yy=y[1,:,1]
        zz=y[1,1,:]

        p1 = contourf(xx,yy,qq)
        title!(p1,"Numerical Solution")
        xlabel!(p1,"x")
        ylabel!(p1,"y")
        zlabel!(p1,"z")
 
        qq=reshape(qe,(Nx,Ny,Nz))
        p2 = contourf(xx,yy,qq)
        title!(p2,"Exact Solution")
        xlabel!(p2,"x")
        ylabel!(p2,"y")
        zlabel!(p2,"z")

        figure_1 = plot(p1,p2,layout = (1,2))
        display(figure_1)
    end

    println("Done") #output

end

#----------------------------------#
# Run the main function
#----------------------------------#
main()
