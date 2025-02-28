#=
---------------------------------------------------------------------
This function plots the solution
Written by F.X. Giraldo on April 19 2019
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function plot_solution(q,coord_CG,periodicity_CG,DG_to_CG,space_method,time,Nx,Ny)

    #Allocate arrays
    Npoin_CG=Nx*Ny
    Npoin=length(q)
    q_sol=zeros(DFloat,Npoin_CG)
    lhowm=zeros(Int64,Npoin_CG)

    #Compute gridpoint solution
    for i=1:Npoin
        ip_CG=periodicity_CG[DG_to_CG[i]]
        q_sol[ip_CG] += q[i]
        lhowm[ip_CG] += 1
    end
    for i=1:Npoin_CG
        q_sol[i]=q_sol[i]/lhowm[i]
    end

    #Plot Solution
    x=reshape(coord_CG[1,:],(Nx,Ny))
    y=reshape(coord_CG[2,:],(Nx,Ny))
    qq=reshape(q_sol,(Nx,Ny))
    xx=x[:,1]
    yy=y[1,:]

    p1 = contourf(xx,yy,qq')
    timec=round(time/(2*π); digits=2)
    title!(p1,"$space_method solution at time = $timec")
    xlabel!(p1,"x")
    ylabel!(p1,"y")
    figure_1=plot(p1)
    display(figure_1)

end
