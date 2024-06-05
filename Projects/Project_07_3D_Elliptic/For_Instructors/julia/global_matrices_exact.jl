#=---------------------------------------------------------------------
This code computes the global matrices using DSS
Written by F.X. Giraldo on April 24, 2019
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216

Cost is O(N^{3d}) for d=3 O(N^9)
---------------------------------------------------------------------=#
function global_matrices_exact(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

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

            #Loop through I points = Rows of Matrices
            for ki=1:Np, ji=1:Np, ii=1:Np
                I=inode[ii,ji,ki]
                Ψ_I=ψ[ii,i]*ψ[ji,j]*ψ[ki,k]
                Ψ_ξ=dψ[ii,i]*ψ[ji,j]*ψ[ki,k]
                Ψ_η=ψ[ii,i]*dψ[ji,j]*ψ[ki,k]
                Ψ_ζ=ψ[ii,i]*ψ[ji,j]*dψ[ki,k]
                ∂Ψ∂x_I=Ψ_ξ*ξ1 + Ψ_η*η1 + Ψ_ζ*ζ1
                ∂Ψ∂y_I=Ψ_ξ*ξ2 + Ψ_η*η2 + Ψ_ζ*ζ2
                ∂Ψ∂z_I=Ψ_ξ*ξ3 + Ψ_η*η3 + Ψ_ζ*ζ3
                
                #Loop through J points = Cols of Matrices
                for kj=1:Np, jj=1:Np, ij=1:Np
                    J=inode[ij,jj,kj]
                    Ψ_J=ψ[ij,i]*ψ[jj,j]*ψ[kj,k]
                    Ψ_ξ=dψ[ij,i]*ψ[jj,j]*ψ[kj,k]
                    Ψ_η=ψ[ij,i]*dψ[jj,j]*ψ[kj,k]
                    Ψ_ζ=ψ[ij,i]*ψ[jj,j]*dψ[kj,k]
                    ∂Ψ∂x_J=Ψ_ξ*ξ1 + Ψ_η*η1 + Ψ_ζ*ζ1
                    ∂Ψ∂y_J=Ψ_ξ*ξ2 + Ψ_η*η2 + Ψ_ζ*ζ2
                    ∂Ψ∂z_J=Ψ_ξ*ξ3 + Ψ_η*η3 + Ψ_ζ*ζ3
                    M[I,J] += wq*Ψ_I*Ψ_J
                    L[I,J] -= wq*( ∂Ψ∂x_I*∂Ψ∂x_J + ∂Ψ∂y_I*∂Ψ∂y_J + ∂Ψ∂z_I*∂Ψ∂z_J )
                end #ij, jj, kj
            end #ii, ji, ki
        end #i, j, k 
    end #e

    return (M,L)
end
