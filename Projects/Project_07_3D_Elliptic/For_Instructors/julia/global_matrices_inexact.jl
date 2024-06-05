#=---------------------------------------------------------------------
This code computes the global matrices using DSS
Written by F.X. Giraldo on April 24, 2019
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216

Updated by F.X. Giraldo on June 4, 2024 to add Inexact Integration

Cost is = N^d[ 9 N^{d-1} ]= 9 N^5 for d=3 instead of O(N^9)
In General we get O( d^d*N^{d+2} )
---------------------------------------------------------------------=#
function global_matrices_inexact(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)
    inode=zeros(Int64,Np,Np,Np)

    @inbounds for e=1:Ne
    
        #Store Coordinates
        for k=1:Np, j=1:Np, i=1:Np
            I=intma[i,j,k,e]
            inode[i,j,k]=I
        end #i,j,k
        
        #Do LGL Integration
        for k=1:Nq, j=1:Nq, i=1:Nq
            wq=ωq[i]*ωq[j]*ωq[k]*jac[i,j,k,e]
            ξ1, ξ2, ξ3 = ξ_x[i,j,k,e], ξ_y[i,j,k,e], ξ_z[i,j,k,e]
            η1, η2, η3 = η_x[i,j,k,e], η_y[i,j,k,e], η_z[i,j,k,e]
            ζ1, ζ2, ζ3 = ζ_x[i,j,k,e], ζ_y[i,j,k,e], ζ_z[i,j,k,e]

            #------------------------Mass Matrix------------------------------#
            I=inode[i,j,k]
            M[I,I] += wq

            #-----------------------Laplacian Matrix--------------------------#            
            #Loop through I points = Rows of Matrix
            for ii=1:Np 
                #ξ derivatives
                I=inode[ii,j,k] #ji=j, ki=k
                Ψ_ξ=dψ[ii,i]*ψ[j,j]*ψ[k,k] 
                ∂Ψ∂x_I=Ψ_ξ*ξ1
                ∂Ψ∂y_I=Ψ_ξ*ξ2
                ∂Ψ∂z_I=Ψ_ξ*ξ3
                
                #Loop through J points = Cols of Matrix
                for ij=1:Np
                    #ξ derivatives
                    J=inode[ij,j,k] #jj=j, kj=k
                    Ψ_ξ=dψ[ij,i]*ψ[j,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_ξ*ξ1
                    ∂Ψ∂y_J=Ψ_ξ*ξ2
                    ∂Ψ∂z_J=Ψ_ξ*ξ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #η derivatives
                    J=inode[i,ij,k] #ii=i, kj=k, swap jj -> ij
                    Ψ_η=ψ[i,i]*dψ[ij,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_η*η1
                    ∂Ψ∂y_J=Ψ_η*η2
                    ∂Ψ∂z_J=Ψ_η*η3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #ζ derivatives
                    J=inode[i,j,ij] #ii=i, jj=j, swap kj -> ij
                    Ψ_ζ=ψ[i,i]*ψ[j,j]*dψ[ij,k]
                    ∂Ψ∂x_J=Ψ_ζ*ζ1
                    ∂Ψ∂y_J=Ψ_ζ*ζ2
                    ∂Ψ∂z_J=Ψ_ζ*ζ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )                    
                end #ij

                #η derivatives
                I=inode[i,ii,k] #ii=i, ki=k, swap ji -> ii
                Ψ_η=ψ[i,i]*dψ[ii,j]*ψ[k,k]
                ∂Ψ∂x_I=Ψ_η*η1
                ∂Ψ∂y_I=Ψ_η*η2
                ∂Ψ∂z_I=Ψ_η*η3
                
                #Loop through J points = Cols of Matrices
                for ij=1:Np
                    #ξ derivatives
                    J=inode[ij,j,k] #jj=j, kj=k
                    Ψ_ξ=dψ[ij,i]*ψ[j,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_ξ*ξ1
                    ∂Ψ∂y_J=Ψ_ξ*ξ2
                    ∂Ψ∂z_J=Ψ_ξ*ξ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #η derivatives
                    J=inode[i,ij,k] #ii=i, kj=k, swap jj -> ij
                    Ψ_η=ψ[i,i]*dψ[ij,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_η*η1
                    ∂Ψ∂y_J=Ψ_η*η2
                    ∂Ψ∂z_J=Ψ_η*η3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #ζ derivatives
                    J=inode[i,j,ij] #ii=i, jj=j, swap kj -> ij
                    Ψ_ζ=ψ[i,i]*ψ[j,j]*dψ[ij,k]
                    ∂Ψ∂x_J=Ψ_ζ*ζ1
                    ∂Ψ∂y_J=Ψ_ζ*ζ2
                    ∂Ψ∂z_J=Ψ_ζ*ζ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )                    
                end #ij

                #ζ derivatives
                I=inode[i,j,ii] #ii=i, ji=j, swap ki -> ii
                Ψ_ζ=ψ[i,i]*ψ[j,j]*dψ[ii,k]
                ∂Ψ∂x_I=Ψ_ζ*ζ1
                ∂Ψ∂y_I=Ψ_ζ*ζ2
                ∂Ψ∂z_I=Ψ_ζ*ζ3
                
                #Loop through J points = Cols of Matrices
                for ij=1:Np
                    #ξ derivatives
                    J=inode[ij,j,k] #jj=j, kj=k
                    Ψ_ξ=dψ[ij,i]*ψ[j,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_ξ*ξ1
                    ∂Ψ∂y_J=Ψ_ξ*ξ2
                    ∂Ψ∂z_J=Ψ_ξ*ξ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #η derivatives
                    J=inode[i,ij,k] #ii=i, kj=k, swap jj -> ij
                    Ψ_η=ψ[i,i]*dψ[ij,j]*ψ[k,k]
                    ∂Ψ∂x_J=Ψ_η*η1
                    ∂Ψ∂y_J=Ψ_η*η2
                    ∂Ψ∂z_J=Ψ_η*η3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )

                    #ζ derivatives
                    J=inode[i,j,ij] #ii=i, jj=j, swap kj -> ij
                    Ψ_ζ=ψ[i,i]*ψ[j,j]*dψ[ij,k]
                    ∂Ψ∂x_J=Ψ_ζ*ζ1
                    ∂Ψ∂y_J=Ψ_ζ*ζ2
                    ∂Ψ∂z_J=Ψ_ζ*ζ3
                    #-------------------Laplacian Matrix------------------#
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )                    
                end #ij

            end #ii
        end #i, j, k 
    end #e

    return (M,L)
end
