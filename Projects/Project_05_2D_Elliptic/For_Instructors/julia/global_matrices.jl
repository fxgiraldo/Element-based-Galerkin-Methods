#---------------------------------------------------------------------#
#This code computes the global matrices using DSS
#Written by F.X. Giraldo on April 24, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function global_matrices(Me,Le,intma,Ne,Np,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)

    #Construct Mass and Differentiation Matrices
    for e=1:Ne
        for l=1:Np, k=1:Np
            J=intma[k,l,e]
            for j=1:Np, i=1:Np
                I=intma[i,j,e]
                M[I,J] += Me[i,j,k,l,e]
                L[I,J] -= Le[i,j,k,l,e]
            end #j,i
        end #k,l
    end #e

    return (M,L)
end
