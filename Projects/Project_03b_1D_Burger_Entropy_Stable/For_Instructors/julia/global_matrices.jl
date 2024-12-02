#---------------------------------------------------------------------#
#This function computes the global Mass and Differentiation Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function global_matrices(Me,De,Dwe,Fe,intma,coord,Npoin,Ne,Np,periodicity,DFloat)

    #Initialize
    M=zeros(DFloat,Npoin,Npoin)
    D=zeros(DFloat,Npoin,Npoin)
    Dw=zeros(DFloat,Npoin,Npoin)
    F=zeros(DFloat,Npoin,Npoin)
    inode=zeros(Int64,Np)
    x=zeros(DFloat,Np)

    #Loop through elements
    for e=1:Ne
        #Store Coordinates
        for i=1:Np
            I=intma[i,e]
            inode[i]=periodicity[I]
            x[i]=coord[I]
        end

        #Element Length & Jacobian
        dx=x[Np]-x[1]
        jac=dx/2

        #Perform DSS
        for j=1:Np
            J=inode[j]
            for i=1:Np
                I=inode[i]
                M[I,J]+=jac*Me[i,j]
                D[I,J]+=De[i,j]
                Dw[I,J]+=Dwe[i,j]
                F[I,J]+=Fe[i,j]
            end #j
        end #i
    end #e

    #Periodicity
    if periodicity[Npoin] == periodicity[1]
        M[Npoin,Npoin]=1
    end

    #check Fmatrix
#     for i=1:Npoin
#         for j=1:Npoin
#             println("i = ",i," j = ",j," F = ",F[i,j])
#         end
#     end

    #Return arrays
    return(M,D,Dw,F)
end
