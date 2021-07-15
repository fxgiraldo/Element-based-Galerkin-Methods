#=
---------------------------------------------------------------------
This function computes the data structures required for a CG or DG method.
Examples include: the coord, intma, and periodicity arrays

Written by F.X. Giraldo on 7/14/2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_CGDG_Storage(space_method,Npoin_CG,Npoin_DG,coord_CG,intma_CG,periodicity_CG,Ne,Np,DFloat)

    if (space_method == "CG")
        Npoin=Npoin_CG
    elseif (space_method == "DG")
        Npoin=Npoin_DG
    end
    DG_to_CG=zeros(Int64,Npoin)

    #Initialize Global Arrays
    coord=zeros(DFloat,2,Npoin)
    intma=zeros(Int64,Np,Np,Ne)
    periodicity=zeros(Int64,Npoin,1)

    if (space_method == "CG")
        coord=copy(coord_CG)
        intma=copy(intma_CG)
        periodicity = copy(periodicity_CG)
        for i=1:Npoin_CG
            DG_to_CG[i]=i
        end
        
    elseif (space_method == "DG")
        Npoin=Npoin_DG
        i_DG=0

        #Get INTMA
        for e=1:Ne
            for j=1:Np, i=1:Np
                i_DG=i_DG + 1
                i_CG=intma_CG[i,j,e]
                intma[i,j,e]=i_DG
                DG_to_CG[i_DG]=i_CG
            end
        end
        
        if (i_DG != Npoin_DG)
            println("ERROR in CREATE_CGDG_STORAGE: i_DG != Npoin_DG = ",i_DG,", ",Npoin_DG)
            exit
        end
        
        #Get COORD
        for i_DG=1:Npoin_DG
            i_CG=DG_to_CG[i_DG]
            coord[1,i_DG]=coord_CG[1,i_CG]
            coord[2,i_DG]=coord_CG[2,i_CG]
        end
        
        #Get periodicity
        for i_DG=1:Npoin_DG
            periodicity[i_DG]=i_DG
        end
        
    end #if

    return(coord,intma,periodicity,DG_to_CG,Npoin)
end

