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

include("global_matrices_exact.jl")
include("global_matrices_inexact.jl")

function global_matrices(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)

    if (Np==Nq)
        (M,L) = global_matrices_inexact(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)
    elseif (Np!=Nq)
        (M,L) = global_matrices_exact(ψ,dψ,ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac,ωq,intma,Ne,Np,Nq,Npoin,DFloat)
    end
    return (M,L)

end
