#---------------------------------------------------------------------#
#This function computes the high-order grid & elements in 3D.
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#

using Plots
include("warp_grid.jl")

function create_grid(Np,Npoin,Ne,Nboun,Nex,Ney,Nez,ξ,plot_grid,warp_grid,DFloat)

    #Initialize Global Arrays
    coord=zeros(DFloat,3,Npoin)
    intma=zeros(Int64,Np,Np,Np,Ne)
    iboun=zeros(Int64,Nboun,1)
 
    #Set some constants
    xmin=-1
    xmax=+1
    ymin=-1
    ymax=+1
    zmin=-1
    zmax=+1
    dx=(xmax-xmin)/Nex
    dy=(ymax-ymin)/Ney
    dz=(zmax-zmin)/Nez
    N=Np-1
    Nx=Nex*N + 1
    Ny=Ney*N + 1
    Nz=Nez*N + 1

    #Initialize Local Arrays
    node=zeros(Int64,Nx,Ny,Nz)

    #GENERATE COORD
    I=0
    kk=0
    for m=1:Nez
        z0=zmin + (m-1)*dz
        if (m == 1)
            n1=1
        else
            n1=2
        end
        for n=n1:Np
            z=( ξ[n]+1 )*dz/2 + z0
            kk=kk+1
            jj=0
            for k=1:Ney
                y0=ymin + (k-1)*dy
                if (k == 1)
                    l1=1
                else
                    l1=2
                end
                for l=l1:Np
                    y=( ξ[l]+1 )*dy/2 + y0
                    jj=jj+1
                    ii=0
                    for i=1:Nex
                        x0=xmin + (i-1)*dx
                        if (i == 1)
                            j1=1
                        else
                            j1=2
                        end
                        for j=j1:Np
                            ii=ii + 1
                            I=I + 1
                            x=( ξ[j]+1 )*dx/2 + x0
                            coord[1,I]=x
                            coord[2,I]=y
                            coord[3,I]=z
                            node[ii,jj,kk]=I
                        end #j
                    end #i
                end #l
            end #k
        end #m 
    end #n

    #GENERATE INTMA
    e=0
    for m=1:Nez
        for k=1:Ney
            for i=1:Nex
                e=e+1
                for n=1:Np
                    kk=N*(m-1) + n
                    for l=1:Np
                        jj=N*(k-1) + l
                        for j=1:Np
                            ii=N*(i-1) + j
                            I=node[ii,jj,kk]
                            intma[j,l,n,e]=I
                        end #j
                    end #l
                end #n
            end #i
        end #k
    end #m

    #GENERATE Boundary Conditions
    B=0

    #Front Face [y=-1]
    j=1
    for k=1:Nz
        for i=1:Nx
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end

    #Back Face [y=+1]
    j=Ny
    for k=1:Nz
        for i=1:Nx
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end

    #Left Face [x=-1]
    i=1
    for k=1:Nz
        for j=2:Ny-1
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end

    #Right Face [x=+1]
    i=Nx
    for k=1:Nz
        for j=2:Ny-1
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end

    #Bottom Face [z=-1]
    k=1
    for j=2:Ny-1
        for i=2:Nx-1
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end

    #Top Face [z=-1]
    k=Nz
    for j=2:Ny-1
        for i=2:Nx-1
            I=node[i,j,k]
            B=B + 1
            iboun[B]=I
        end
    end
   
    if (B != Nboun)
        disp("Error in Create_Grid")
        return
    end

    #Warp Grid
    if (warp_grid)
        warp_grid!(coord,Npoin,DFloat)
    end

    #Plot Grid
    title_text=["Grid : Ne =  $Ne,  N = $N"]
    if (plot_grid)
        x=zeros(DFloat,5,1)
        y=zeros(DFloat,5,1)
        z=zeros(DFloat,5,1)
        x[1]=-1; y[1]=-1; z[1]=-1
        x[2]=+1; y[2]=-1; z[2]=-1
        x[3]=+1; y[3]=+1; z[3]=-1
        x[4]=-1; y[4]=+1; z[4]=-1
        x[5]=-1; y[5]=-1; z[5]=-1
        p1=plot(x,y,z,linecolor="black",legend=false)

        for e=1:Ne

            #Build x-y plane grid for each z value
            for k=1:Np
                for j=1:N
                    for i=1:N
                        i1=intma[i,j,k,e]
                        i2=intma[i+1,j,k,e]
                        i3=intma[i+1,j+1,k,e]
                        i4=intma[i,j+1,k,e]
                        x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
                        x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
                        x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
                        x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
                        x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
                        p1=plot!(p1,x,y,z,linecolor="red",legend=false)
                    end
                end
            end

            #Build x-z plane grid for each y value
            for j=1:Np
                for k=1:N
                    for i=1:N
                        i1=intma[i,j,k,e]
                        i2=intma[i+1,j,k,e]
                        i3=intma[i+1,j,k+1,e]
                        i4=intma[i,j,k+1,e]
                        x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
                        x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
                        x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
                        x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
                        x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
                        p1=plot!(p1,x,y,z,linecolor="red",legend=false)
                    end
                end
            end

            #Build y-z plane grid for each x value
            for i=1:Np
                for k=1:N
                    for j=1:N
                        i1=intma[i,j,k,e]
                        i2=intma[i,j+1,k,e]
                        i3=intma[i,j+1,k+1,e]
                        i4=intma[i,j,k+1,e]
                        x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
                        x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
                        x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
                        x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
                        x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
                        p1=plot!(p1,x,y,z,linecolor="red",legend=false)
                   end
                end
            end

            #---------------------Build Element Boundary----------------------#
            #Build x-y plane grid for k=1 value
            i1=intma[1,1,1,e]
            i2=intma[Np,1,1,e]
            i3=intma[Np,Np,1,e]
            i4=intma[1,Np,1,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)

            #Build x-y plane grid for k=Np value
            i1=intma[1,1,Np,e]
            i2=intma[Np,1,Np,e]
            i3=intma[Np,Np,Np,e]
            i4=intma[1,Np,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)

            #Build x-z plane grid for j=1
            i1=intma[1,1,1,e]
            i2=intma[Np,1,1,e]
            i3=intma[Np,1,Np,e]
            i4=intma[1,1,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)

            #Build x-z plane grid for j=Np
            i1=intma[1,Np,1,e]
            i2=intma[Np,Np,1,e]
            i3=intma[Np,Np,Np,e]
            i4=intma[1,Np,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)

            #Build y-z plane grid for each i=1 value
            i1=intma[1,1,1,e]
            i2=intma[1,Np,1,e]
            i3=intma[1,Np,Np,e]
            i4=intma[1,1,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)

            #Build y-z plane grid for each i=Np value
            i1=intma[Np,1,1,e]
            i2=intma[Np,Np,1,e]
            i3=intma[Np,Np,Np,e]
            i4=intma[Np,1,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]; z[1]=coord[3,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]; z[2]=coord[3,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]; z[3]=coord[3,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]; z[4]=coord[3,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]; z[5]=coord[3,i1]
            p1=plot!(p1,x,y,z,linecolor="blue",legend=false)
        end #e
        display(plot(p1))
    end #plot_grid

    return(coord,intma,iboun)
end
