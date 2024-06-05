#=
---------------------------------------------------------------------
This function computes the element mass and Laplacian matrices.

The strategy follows Algorithm 12.7 in F.X. Giraldo's Introduction to Element-based 
Galerkin Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.

It combines the Elemental and Global into one Algorithm

Written by F.X. Giraldo on 5/2024
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216

Cost is O(N^{3d})           
---------------------------------------------------------------------
=#

function global_matrices_exact(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)

    #Construct Mass and Differentiation Matrices
    for e=1:Ne
        #Quadrature Points
        for j=1:Nq, i=1:Nq
            wq=ωq[i]*ωq[j]*jac[i,j,e]
            e_x=ξ_x[i,j,e]; e_y=ξ_y[i,j,e]
            n_x=η_x[i,j,e]; n_y=η_y[i,j,e]
            
            #I points = Rows of the Matrix
            for ji=1:Np, ii=1:Np
                I=intma[ii,ji,e]
                Ψ_I=ψ[ii,i]*ψ[ji,j]
                Ψ_ξ=dψ[ii,i]*ψ[ji,j]
                Ψ_η=ψ[ii,i]*dψ[ji,j]
                dΨdx_I=Ψ_ξ*e_x + Ψ_η*n_x
                dΨdy_I=Ψ_ξ*e_y + Ψ_η*n_y

                #J points = Cols of the Matrix
                for jj=1:Np, ij=1:Np
                    J=intma[ij,jj,e]
                    Ψ_J=ψ[ij,i]*ψ[jj,j]
                    Ψ_ξ=dψ[ij,i]*ψ[jj,j]
                    Ψ_η=ψ[ij,i]*dψ[jj,j]
                    dΨdx_J=Ψ_ξ*e_x + Ψ_η*n_x
                    dΨdy_J=Ψ_ξ*e_y + Ψ_η*n_y
                    #--------------------Mass and Laplacian Matrices----------------#
                    M[I,J]+=wq*Ψ_I*Ψ_J
                    L[I,J]-=wq*( dΨdx_I*dΨdx_J + dΨdy_I*dΨdy_J )
                end #ij,jj
            end #ii,ji
        end #i,j
    end #e

    return (M,L)
end
