#---------------------------------------------------------------------#
#This function computes the high-order grid & elements in 2D.
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#

using Plots
include("warp_grid.jl")
include("rotate_grid.jl")

function create_grid(Np,Npoin,Ne,Nboun,Nex,Ney,ξ,plot_grid,warp_grid,rotate_grid,DFloat)

    #Initialize Global Arrays
    coord=zeros(DFloat,2,Npoin)
    intma=zeros(Int64,Np,Np,Ne)
    bsido=zeros(Int64,4,Nboun)
    periodicity=zeros(Int64,Npoin,1)

    #Initialize Local Arrays
    node=zeros(Int64,Npoin,Npoin)

    #Set some constants
    xmin=-1
    xmax=+1
    ymin=-1
    ymax=+1
    dx=(xmax-xmin)/Nex
    dy=(ymax-ymin)/Ney
    N=Np-1
    nx=Nex*N + 1
    ny=Ney*N + 1

    #GENERATE COORD
    I=0
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
                    node[ii,jj]=I
                end #j

            end #i
        end #l
    end #k

    #GENERATE INTMA
    e=0
    for k=1:Ney
        for i=1:Nex
            e=e+1
            for l=1:Np
                jj=N*(k-1) + l
                for j=1:Np
                    ii=N*(i-1) + j
                    I=node[ii,jj]
                    intma[j,l,e]=I
                end #j
            end #l
        end #i
    end #k

    #Generate BSIDO
    b=0
    #Bottom Boundary
    for i=1:Nex
        e=i
        b=b+1
        i1=(i-1)*(Np-1) + 1
        i2=(i-1)*(Np-1) + Np
        I1=node[i1,1]
        I2=node[i2,1]
        bsido[1,b]=I1
        bsido[2,b]=I2
        bsido[3,b]=e
        bsido[4,b]=6
    end

    #Right Boundary
    for i=1:Ney
        e=(Nex)*(i)
        b=b+1
        i1=(i-1)*(Np-1) + 1
        i2=(i-1)*(Np-1) + Np
        I1=node[nx,i1]
        I2=node[nx,i2]
        bsido[1,b]=I1
        bsido[2,b]=I2
        bsido[3,b]=e
        bsido[4,b]=6
    end 

    #Top Boundary
    for i=Nex:-1:1 
        e=Ne - (Nex - i)
        b=b+1
        i1=(i-1)*(Np-1) + Np
        i2=(i-1)*(Np-1) + 1
        I1=node[i1,ny]
        I2=node[i2,ny]
        bsido[1,b]=I1
        bsido[2,b]=I2
        bsido[3,b]=e
        bsido[4,b]=6
    end 

    #Left Boundary
    for i=Ney:-1:1
        e=(Nex)*(i-1) + 1
        b=b+1
        i1=(i-1)*(Np-1) + Np
        i2=(i-1)*(Np-1) + 1
        I1=node[1,i1]
        I2=node[1,i2]
        bsido[1,b]=I1
        bsido[2,b]=I2
        bsido[3,b]=e
        bsido[4,b]=6
    end

    if (b != Nboun)
        println("Error in Create_Grid: b != Nboun =>  Nboun = ",Nboun," b = ",b)
        exit
    end

    #Periodicity
    for i=1:Npoin
        periodicity[i]=i;
    end

    #X-Periodicity
    for i=1:ny
        i1=node[1,i]
        i2=node[nx,i]
        periodicity[i2]=i1
    end

    #Y-Periodicity
    for i=1:nx
        i1=node[i,1]
        i2=node[i,ny]
        periodicity[i2]=periodicity[i1]
    end

    #Warp Grid
    if (warp_grid)
        warp_grid!(coord,Npoin,DFloat)
    end

    #Rotate Grid
    if (rotate_grid)
        rotate_grid!(coord,Npoin,DFloat)
    end

    #Plot Grid
    title_text=["Grid : Ne =  $Ne,  N = $N"]
    if (plot_grid)
        x=zeros(DFloat,5,1)
        y=zeros(DFloat,5,1)
        x[1]=-1; y[1]=-1
        x[2]=+1; y[2]=-1
        x[3]=+1; y[3]=+1
        x[4]=-1; y[4]=+1
        x[5]=-1; y[5]=-1
        p1=plot(x,y,linecolor="black",legend=false)

        for e=1:Ne
            for j=1:N
                for i=1:N
                    i1=intma[i,j,e]
                    i2=intma[i+1,j,e]
                    i3=intma[i+1,j+1,e]
                    i4=intma[i,j+1,e]
                    x[1]=coord[1,i1]; y[1]=coord[2,i1]
                    x[2]=coord[1,i2]; y[2]=coord[2,i2]
                    x[3]=coord[1,i3]; y[3]=coord[2,i3]
                    x[4]=coord[1,i4]; y[4]=coord[2,i4]
                    x[5]=coord[1,i1]; y[5]=coord[2,i1]
                    p1=plot!(p1,x,y,linecolor="red",legend=false)
                end
            end
            i1=intma[1,1,e]
            i2=intma[Np,1,e]
            i3=intma[Np,Np,e]
            i4=intma[1,Np,e]
            x[1]=coord[1,i1]; y[1]=coord[2,i1]
            x[2]=coord[1,i2]; y[2]=coord[2,i2]
            x[3]=coord[1,i3]; y[3]=coord[2,i3]
            x[4]=coord[1,i4]; y[4]=coord[2,i4]
            x[5]=coord[1,i1]; y[5]=coord[2,i1]
            p1=plot!(p1,x,y,linecolor="blue",legend=false)
        end
        display(plot(p1))
    end

    return(coord,intma,bsido,periodicity)
end
