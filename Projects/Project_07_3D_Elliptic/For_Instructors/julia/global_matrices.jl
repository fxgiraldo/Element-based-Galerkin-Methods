#---------------------------------------------------------------------#
#This code computes the global matrices using DSS
#Written by F.X. Giraldo on April 24, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function global_matrices(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)
    inode=zeros(Int64,Np,Np,Np)

    for e=1:Ne
    
        #Store Coordinates
        for k=1:Np, j=1:Np, i=1:Np
            I=intma[i,j,k,e]
            inode[i,j,k]=I
        end #i,j,k
        
        #Do LGL Integration
        for n=1:Nq, m=1:Nq, l=1:Nq
            wq=ωq[l]*ωq[m]*ωq[n]*jac[l,m,n,e]
            
            #Loop through I points
            for k=1:Np, j=1:Np, i=1:Np
                I=inode[i,j,k]
                h_i=ψ[i,l]*ψ[j,m]*ψ[k,n]
                h_ξ=dψ[i,l]*ψ[j,m]*ψ[k,n]
                h_η=ψ[i,l]*dψ[j,m]*ψ[k,n]
                h_ζ=ψ[i,l]*ψ[j,m]*dψ[k,n]
                δhδx_i=h_ξ*ξ_x[l,m,n,e] + h_η*η_x[l,m,n,e] + h_ζ*ζ_x[l,m,n,e]
                δhδy_i=h_ξ*ξ_y[l,m,n,e] + h_η*η_y[l,m,n,e] + h_ζ*ζ_y[l,m,n,e]
                δhδz_i=h_ξ*ξ_z[l,m,n,e] + h_η*η_z[l,m,n,e] + h_ζ*ζ_z[l,m,n,e]
                
                #Loop through J points
                for kj=1:Np, jj=1:Np, ij=1:Np
                    J=inode[ij,jj,kj]
                    h_j=ψ[ij,l]*ψ[jj,m]*ψ[kj,n]
                    h_ξ=dψ[ij,l]*ψ[jj,m]*ψ[kj,n]
                    h_η=ψ[ij,l]*dψ[jj,m]*ψ[kj,n]
                    h_ζ=ψ[ij,l]*ψ[jj,m]*dψ[kj,n]
                    δhδx_j=h_ξ*ξ_x[l,m,n,e] + h_η*η_x[l,m,n,e] + h_ζ*ζ_x[l,m,n,e]
                    δhδy_j=h_ξ*ξ_y[l,m,n,e] + h_η*η_y[l,m,n,e] + h_ζ*ζ_y[l,m,n,e]
                    δhδz_j=h_ξ*ξ_z[l,m,n,e] + h_η*η_z[l,m,n,e] + h_ζ*ζ_z[l,m,n,e]
                    M[I,J] += wq*h_i*h_j
                    L[I,J] -= wq*( δhδx_i*δhδx_j + δhδy_i*δhδy_j + δhδz_i*δhδz_j )
                end #ij, jj, kj
            end #i, j, k
        end #l, m, n 
    end #e

    return (M,L)
end
