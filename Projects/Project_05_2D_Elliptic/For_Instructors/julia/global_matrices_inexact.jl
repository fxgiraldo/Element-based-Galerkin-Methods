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

Cost is = N^d[ 4N^{d-1} ]= 2N^4 for d=2

In General we get O( d^d*N^{d+2} )
---------------------------------------------------------------------
=#

function global_matrices_inexact(ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    #Initialize Matrices
    M=zeros(DFloat,Npoin,Npoin)
    L=zeros(DFloat,Npoin,Npoin)

    #Construct Mass and Differentiation Matrices
    for e=1:Ne
        for j=1:Nq, i=1:Nq
            wq=ωq[i]*ωq[j]*jac[i,j,e]
            e_x=ξ_x[i,j,e]; e_y=ξ_y[i,j,e]
            n_x=η_x[i,j,e]; n_y=η_y[i,j,e]

            #------------------------Mass Matrix------------------------------#
            I=intma[i,j,e]
            M[I,I]+=wq

            #-----------------------Laplacian Matrix--------------------------#            
            #Loop through I points = Rows of the Matrix
            for ii=1:Np
     
                #XI derivatives
                I=intma[ii,j,e] #ji=j
                Ψ_ξ=dψ[ii,i]*ψ[j,j] #ji=j
                dΨdx_I=Ψ_ξ*e_x
                dΨdy_I=Ψ_ξ*e_y

                #Loop through J points = Cols of the Matrix
                for ij=1:Np
                    J=intma[ij,j,e] #jj=j
                    Ψ_ξ=dψ[ij,i]*ψ[j,j] #jj=j
                    dΨdx_J=Ψ_ξ*e_x
                    dΨdy_J=Ψ_ξ*e_y
                    L[I,J]=L[I,J] - wq*(dΨdx_I*dΨdx_J + dΨdy_I*dΨdy_J);
                    J=intma[i,ij,e] #ii=i and swap jj-> ij
                    Ψ_η=ψ[i,i]*dψ[ij,j] #ii=i and swap jj-> ij
                    dΨdx_J=Ψ_η*n_x
                    dΨdy_J=Ψ_η*n_y
                    L[I,J]=L[I,J] - wq*(dΨdx_I*dΨdx_J + dΨdy_I*dΨdy_J);
                end #ij
     
                #ETA derivatives
                I=intma[i,ii,e] #ii=i and swapped ji->ii
                Ψ_η=ψ[i,i]*dψ[ii,j] #ii=i and swapped ji->ii
                dΨdx_I=Ψ_η*n_x
                dΨdy_I=Ψ_η*n_y

                #Loop through J points = Cols of the Matrix
                for ij=1:Np
                    J=intma[ij,j,e] #jj=j
                    Ψ_ξ=dψ[ij,i]*ψ[j,j] #jj=j
                    dΨdx_J=Ψ_ξ*e_x
                    dΨdy_J=Ψ_ξ*e_y
                    L[I,J]=L[I,J] - wq*(dΨdx_I*dΨdx_J + dΨdy_I*dΨdy_J);
                    J=intma[i,ij,e] #ii=i and swap jj-> ij
                    Ψ_η=ψ[i,i]*dψ[ij,j] #ii=i and swap jj-> ij
                    dΨdx_J=Ψ_η*n_x
                    dΨdy_J=Ψ_η*n_y
                    L[I,J]=L[I,J] - wq*(dΨdx_I*dΨdx_J + dΨdy_I*dΨdy_J);
                end #ij
            end #ii
     
        end #i,j
    end #e

    return (M,L)
end
