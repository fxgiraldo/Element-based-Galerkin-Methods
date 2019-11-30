#---------------------------------------------------------------------#
#This code computes the element mass and differentiation matrices
#Written by F.X. Giraldo on April 24, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function element_matrices(ψ,dψ,Np,Nq,ωq,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Np,Np)
    D=zeros(DFloat,Np,Np)
    Dw=zeros(DFloat,Np,Np)
    F=zeros(DFloat,Np,Np)

    #Construct Mass and Differentiation Matrices
    @inbounds for k=1:Nq, j=1:Np, i=1:Np
        M[i,j]+=ωq[k]*ψ[i,k]*ψ[j,k]
        D[i,j]+=ωq[k]*ψ[i,k]*dψ[j,k]
        Dw[i,j]+=ωq[k]*dψ[i,k]*ψ[j,k]
    end
    F[1,1]=-1; F[Np,Np]=+1

    return (M,D,Dw,F)
end
