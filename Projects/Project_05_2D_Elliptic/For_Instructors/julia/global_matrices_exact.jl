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
---------------------------------------------------------------------
=#

function global_matrices_exact(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)

    #Construct Mass and Differentiation Matrices
    for e=1:Ne
        for l=1:Nq, k=1:Nq
            wq=ωq[k]*ωq[l]*jac[k,l,e]
            for j=1:Np, i=1:Np
                I=intma[i,j,e]
                Ψ_JK=ψ[i,k]*ψ[j,l] #h_ik*h_jl
                dΨdx_JK=dψ[i,k]*ψ[j,l]*ξ_x[k,l,e] + ψ[i,k]*dψ[j,l]*η_x[k,l,e]
                dΨdy_JK=dψ[i,k]*ψ[j,l]*ξ_y[k,l,e] + ψ[i,k]*dψ[j,l]*η_y[k,l,e]
                for n=1:Np, m=1:Np
                    J=intma[m,n,e]
                    Ψ_IK=ψ[m,k]*ψ[n,l] #h_ik*h_jl
                    dΨdx_IK=dψ[m,k]*ψ[n,l]*ξ_x[k,l,e] + ψ[m,k]*dψ[n,l]*η_x[k,l,e]
                    dΨdy_IK=dψ[m,k]*ψ[n,l]*ξ_y[k,l,e] + ψ[m,k]*dψ[n,l]*η_y[k,l,e]
                    M[I,J]+=wq*Ψ_IK*Ψ_JK
                    L[I,J]-=wq*( dΨdx_IK*dΨdx_JK + dΨdy_IK*dΨdy_JK )
                end #m,n
            end #j,i
        end #k,l
    end #e

    return (M,L)
end
