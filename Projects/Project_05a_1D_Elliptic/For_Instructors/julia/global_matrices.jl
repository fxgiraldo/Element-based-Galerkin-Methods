#---------------------------------------------------------------------#
#This function computes the global Mass and Differentiation Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function global_matrices(ψ,dψdξ,ωq,intma,coord,Npoin,Ne,Np,Nq,DFloat)

    #Initialize
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)
    inode=zeros(Int64,Np)
    x=zeros(DFloat,Np)

    #Loop through elements
    for e=1:Ne
        #Store Coordinates
        for i=1:Np
            I=intma[i,e]
            inode[i]=I
            x[i]=coord[I]
        end

        #Element Length & Jacobian
        dx=x[Np]-x[1]
        dxdξ=dx/2
        dξdx=2.0/dx

        #Build Matrices
        @inbounds for k=1:Nq
            for j=1:Np
                J=inode[j]
                for i=1:Np
                    I=inode[i]
                    M[I,J]+=ωq[k]*ψ[i,k]*ψ[j,k]*dxdξ
                    L[I,J]+=ωq[k]*(dψdξ[i,k]*dξdx)*(dψdξ[j,k]*dξdx)*dxdξ
                 end #i
            end #j
        end #k
    end #e

    #Impose Dirichlet BC on L
    for J=1:Npoin
        L[1,J]=0
        L[Npoin,J]=0
    end
    L[1,1]=DFloat(1.0)
    L[Npoin,Npoin]=DFloat(1.0)
    
    #Return arrays
    return(M,L)
end
