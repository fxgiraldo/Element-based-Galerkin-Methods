#=
-------------------------------------------------------------------------------------------------------------
This file runs the 2D Wave Equation using the unified CG/DG method with AGGP storage 
introduced in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.

It uses the LSRK45 method as the time-integrator and Rusanov fluxes for DG.

Written by F.X. Giraldo on July 10, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943

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
using LinearAlgebra

include("QuadraturePoints.jl")
include("element_matrices.jl")
include("global_matrices.jl")
include("create_grid.jl")
include("create_CGDG_storage.jl")
include("create_side.jl")
include("create_face.jl")
include("compute_normals.jl")
include("create_face_periodicity.jl")
include("compute_metrics.jl")
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norm.jl")
include("ti_LSRK.jl")

#Some Constants
DFloat = Float64
machine_zero=eps(DFloat)

#Main Function
function main()

    #------------------------Only Change stuff here-------------------#
    N=8 #polynomial order
    Nex=4 #number of elements along x (ξ)
    Ney=4 #number of elements along y (η)
    space_method="DG" #CG or DG
    time_final=DFloat(0.25) #in revolutions
    plot_grid=false
    plot_solution=true
    plot_movie=false
    warp_grid=false
    #------------------------Only Change stuff here-------------------#

    Q=N #changing this requires modifying CREATE_RHS_FLUX which assumes Q=N
    Np=N+1
    Nq=Q+1
    Npts=Np^2
    Nqs=Nq^2
    ipoints=1 #lobatto
    qpoints=1 #lobatto
    Courant_max=DFloat(0.5)
    icase=1 #CCW moving Gaussian
    Ne=Nex*Ney
    Nx=Nex*N+1
    Ny=Ney*N+1
    Npoin_CG=Nx*Ny
    Npoin_DG=Npts*Ne
    Npoin=Npoin_CG
    Nboun=2*Nex + 2*Ney
    Nface=2*Ne + Nex + Ney
    if icase == 1
        c=2*π;
    end
    rotate_grid=false #not working with periodicity
    time_final=time_final*c
    @show(N,Q,Ne,Npoin_CG,Npoin_DG,Nboun,Nface,Nx,Ny)

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
    (coord_CG,intma_CG,bsido_CG,periodicity_CG) = create_grid(Np,Npoin_CG,Ne,Nboun,Nex,Ney,ξ,plot_grid,warp_grid,rotate_grid,DFloat)

    #Store the CGDG-Storage Grid
    (coord,intma,periodicity,DG_to_CG,Npoin) = create_CGDG_Storage(space_method,Npoin_CG,Npoin_DG,coord_CG,intma_CG,periodicity_CG,Ne,Np,DFloat)
    @show(space_method,Npoin)

    #Construct Metric Terms
    (ξ_x,ξ_y,η_x,η_y,jac) = compute_metrics(coord_CG,intma_CG,ψ,dψ,Ne,Np,Nq,DFloat)

    #Construct Face arrays
    (iside,~) = create_side(intma_CG,bsido_CG,Npoin_CG,Ne,Nboun,Nface,Np,DFloat)
    (face,mapL,mapR) = create_face(iside,intma_CG,Nface,Np)
    (normals,jac_face) = compute_normals(face,intma_CG,coord_CG,Nface,Np,Nq,ψ,dψ)
    (face,iside) = create_face_periodicity!(face,iside,coord_CG,Nface,Nboun)

    #Construct Element Matrices
    (Me,De_x,De_y) = element_matrices(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,Ne,Np,Nq,DFloat)

    #Construct Global Matrices
    M = global_matrices(Me,intma,periodicity,Ne,Np,Npoin,DFloat)

    #Compute Exact Solution
    time=DFloat(0.0)
    (qe,ue) = exact_solution(coord,Npoin,time,icase,DFloat)
    q0=copy(qe)

    #Recompute DT, NTIME
    dx_min=coord[1,2]-coord[1,1]
    umax=maximum(ue)
    dt=Courant_max*dx_min/umax
    courant=umax*dt/dx_min
    ntime = ceil(Int64, time_final / dt)
    println("Courant = ",courant," dt = ",dt," ntime = ",ntime)

    #--------------------------------------------------------------------------------------------------------
    #Students Add their Contributions in CREATE_RHS which is called inside of ti_LSRK
    #--------------------------------------------------------------------------------------------------------

    #Time Integration
    (q0,time) = ti_LSRK!(q0,ue,coord,M,De_x,De_y,intma,periodicity,face,mapL,mapR,normals,jac_face,ωq,time,ntime,dt,space_method,plot_movie,Nx,Ny,coord_CG,periodicity_CG,DG_to_CG,DFloat)
    println(" itime = ",ntime," time = ",time," qmax = ",maximum(q0)," qmin = ",minimum(q0))

    #Compute L2 Norm
    (qe,ue) = exact_solution(coord,Npoin,time,icase,DFloat)
    norm2 = compute_norm(q0,qe,Npoin,DFloat)
    println("norm.L2 = ",norm2)

    #Plot Solution
    if (plot_solution)
        #Compute gridpoint solution
        q_sol=zeros(DFloat,Npoin_CG)
        qe_sol=zeros(DFloat,Npoin_CG)
        lhowm=zeros(Int64,Npoin_CG)
        for i=1:Npoin
            ip_CG=periodicity_CG[DG_to_CG[i]]
            q_sol[ip_CG] += q0[i]
            qe_sol[ip_CG] += qe[i]
            lhowm[ip_CG] += 1
        end
        for i=1:Npoin_CG
            q_sol[i]=q_sol[i]/lhowm[i]
            qe_sol[i]=qe_sol[i]/lhowm[i];
        end

        x=reshape(coord_CG[1,:],(Nx,Ny))
        y=reshape(coord_CG[2,:],(Nx,Ny))
        xx=x[:,1]
        yy=y[1,:]

        #Numerical Solution
        qq=reshape(q_sol,(Nx,Ny))
        p1 = contourf(xx,yy,qq')
        title!(p1,"$space_method solution")
        xlabel!(p1,"x")
        ylabel!(p1,"y")

        #Exact Solution
        qq=reshape(qe_sol,(Nx,Ny))
        p2 = contourf(xx,yy,qq')
        title!(p2,"Exact solution")
        xlabel!(p2,"x")
        ylabel!(p2,"y")
        figure_1=plot(p1,p2)
        display(figure_1)
    end

    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
