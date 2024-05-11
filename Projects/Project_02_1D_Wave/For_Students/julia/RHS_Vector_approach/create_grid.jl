#---------------------------------------------------------------------#
#This function computes the Grid and Elements in 1D.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_grid(Np,Npoin_cg,Npoin_dg,Ne,x,DFloat)

    coord_cg=zeros(DFloat,Npoin_cg)
    coord_dg=zeros(DFloat,Npoin_dg)
    intma_cg=zeros(Int64,Np,Ne)
    intma_dg=zeros(Int64,Np,Ne)
    periodicity_cg=zeros(Int64,Npoin_cg)
    periodicity_dg=zeros(Int64,Npoin_dg)

    #Set some constants
    xmin=-1
    xmax=+1
    dx=(xmax-xmin)/Ne

    #Generate Grid Points
    I=1
    coord_cg[1]=xmin
    for e=1:Ne
        x0=xmin + (e-1)*dx
        intma_cg[1,e]=I
        for i=2:Np
            I+=1
            coord_cg[I]=( x[i]+1 )*dx/2 + x0
            intma_cg[i,e]=I
        end
    end

    #Generate Periodicity_CG array
    for I=1:Npoin_cg
        periodicity_cg[I]=I
    end
    periodicity_cg[Npoin_cg]=1

    #Generate INTMA_DG
    I=0
    for e=1:Ne, i=1:Np
        I+=1
        intma_dg[i,e]=I
    end

    #Generate COORD_DG
    for e=1:Ne, i=1:Np
        Icg=intma_cg[i,e]
        Idg=intma_dg[i,e]
        coord_dg[Idg]=coord_cg[Icg]
    end

    #Generate Periodicity_DG array
    for I=1:Npoin_dg
        periodicity_dg[I]=I
    end

    return (coord_cg, coord_dg, intma_cg, intma_dg, periodicity_cg, periodicity_dg)
end
