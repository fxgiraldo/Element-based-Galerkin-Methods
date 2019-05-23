using Plots

include("QuadraturePoints.jl")
include("element_matrices.jl")
include("global_matrices.jl")
include("create_grid.jl")
include("compute_metrics.jl")
include("global_matrices.jl")
include("exact_solution.jl")
include("compute_norm.jl")

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
    N=16
    Q=N+1
    Np=N+1
    Nq=Q+1
    Npts=Np^2
    Nqs=Nq^2
    ipoints=1
    qpoints=1
    space_method="CG"
    c=DFloat(1) #Constant in Exact solution
    Nex=1
    Ney=1
    Ne=Nex*Ney
    Nx=Nex*N+1
    Ny=Ney*N+1
    Npoin=Nx*Ny
    Nboun=2*Nx + 2*(Ny-2);
    plot_grid=true
    plot_solution=false
    warp_grid=false
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
    (coord,intma,iboun) = create_grid(Np,Npoin,Ne,Nboun,Nex,Ney,ξ,plot_grid,warp_grid,DFloat)

    #Construct Metric Terms
    (ξ_x,ξ_y,η_x,η_y,jac) = compute_metrics(coord,intma,ψ,dψ,ωq,Ne,Np,Nq,DFloat)

    #Construct Element Matrices
#    (Me,Le) = element_matrices(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,Ne,Np,Nq,DFloat)

    #Construct Global Matrices
#    (M,L) = global_matrices(Me,Le,intma,Ne,Np,Npoin,DFloat)

    #Compute Exact Solution
    (qe,fe) = exact_solution(coord,Npoin,c,DFloat)

    #Apply Dirichlet BCs: Impose exact solution
    #=
    for i=1:Nboun
        I=iboun[i]
        L[I,:] = ...
        L[I,I] = ...
        R[I]   = ...
    end#i

    #Solve Matrix problem
    q=-L\R

    #Compute L2 Norm
    norm = compute_norm(q,qe,Npoin,DFloat)
    println("2-norm = ",norm[1])
    =#

    #Plot Solution
    if (plot_solution)
        x=reshape(coord[1,:],(Nx,Ny))
        y=reshape(coord[2,:],(Nx,Ny))
        qq=reshape(q,(Nx,Ny))
        xx=x[:,1]
        yy=y[1,:]
        p1=contour(xx,yy,qq,xlabel="x",ylabel="y",zlabel="q")
        p2=contourf(xx,yy,qq,xlabel="x",ylabel="y",zlabel="q")
        p3=surface(xx,yy,qq,xlabel="x",ylabel="y",zlabel="q")
        plot_handle=plot(p1,p2,layout=(2,1),legend=false)
        display(p2)
    end

    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()
